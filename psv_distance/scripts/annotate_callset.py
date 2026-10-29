#!/usr/bin/env python3
"""Annotate official benchmark-side observations after context-aware matching.

Alleles are normalized only in a distance-analysis derivative; the benchmark
observations, BVT classes, decisions and count multiplicities are preserved.
Equal-length substitution alleles use only truly changed reference bases;
unchanged internal MNP bases never create artificial zero-distance overlaps.
The same changed-base rule applies to equal-length PSV alleles. PSV start/end
fields retain source bounding spans; nearest_psv_effective_spans0 records the
actual geometric runs. Nonoverlap sensitivity excludes whole parent PSV loci.
"""
import argparse,csv,gzip,json,os,random,shlex,subprocess,time
from collections import Counter,defaultdict
from pathlib import Path
import pysam
from benchmark_counts import observations,record_regions,passes_filter,compare_counts,COUNT_COLUMNS
from intervals import BedIndex,read_bed
from psv_distance import PSVIndex,allele_effective_spans,span_distance

SCOPES=('Segdup','Segdup_lowmap')
NORMALIZER='/path/to/hap.py/bin/bcftools'

def selected_alleles(record,side):
    gt=record.samples[side].get('GT') or ()
    indices=sorted({i for i in gt if i is not None and i>0 and i<=len(record.alts or ())})
    fallback=not indices
    # A missing/nonvariant GT cannot identify which side carries an ALT.
    # Keep the benchmark observation but label distance unsupported; never borrow the other side's allele.
    return [(i,record.alts[i-1]) for i in indices if record.alts[i-1] not in ('.','*','<NON_REF>','<*>') and all(b in 'ACGTNacgtn' for b in record.alts[i-1])],fallback

def write_counts(path,counter,callset,vtype):
    columns=['callset','tech','caller','scope','policy','vtype','distance_bin','family']+list(COUNT_COLUMNS)+['QUERY.FILTERED']
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader()
        for (scope,policy,bin,family),counts in sorted(counter.items()):
            d={k:callset[k] for k in ('tech','caller')};d.update(callset=callset['id'],scope=scope,policy=policy,vtype=vtype,distance_bin=bin,family=family)
            d.update({c:counts[c] for c in columns[8:]});w.writerow(d)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',default='config/inputs.json');ap.add_argument('--index',type=int);ap.add_argument('--callset')
    ap.add_argument('--vtype',required=True,choices=['SNP','INDEL']);ap.add_argument('--prepared',default='prepared');ap.add_argument('--quantified',default='work/quantified');ap.add_argument('--out-dir',default='results');ap.add_argument('--manual-checks',type=int,default=40)
    a=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):ap.error('Run full-file annotation in a Slurm allocation')
    start=time.time();cfg=json.load(open(a.config))
    callset=cfg['callsets'][a.index] if a.index is not None else next(c for c in cfg['callsets'] if c['id']==a.callset)
    ident=callset['id'];out=Path(a.out_dir);out.mkdir(exist_ok=True,parents=True)
    prefix=out/f'{ident}_{a.vtype}'
    if Path(str(prefix)+'.validation.json').exists():ap.error('Completed output exists; preserve it and use a new out-dir')
    work=Path('work/annotation')/f'{ident}_{a.vtype}';work.mkdir(exist_ok=True,parents=True)
    qprefix=Path(a.quantified)/ident;vcf=Path(str(qprefix)+'.vcf.gz')
    meta=work/'observations.jsonl.gz';allelevcf=work/'alleles.vcf';counts=Counter();qc=Counter();input_allele_ids=set()
    with pysam.FastaFile(cfg['reference']) as fasta,gzip.open(meta,'wt') as mf,open(allelevcf,'w') as af:
        af.write('##fileformat=VCFv4.2\n##reference='+cfg['reference']+'\n')
        for chrom,length in zip(fasta.references,fasta.lengths):af.write(f'##contig=<ID={chrom},length={length}>\n')
        af.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        with pysam.VariantFile(str(vcf)) as source:
            for n,record in enumerate(source):
                qc['benchmark_records_scanned']+=1
                scopes=record_regions(record)&set(SCOPES)
                if not scopes:continue
                obsall=observations(record,'ALL');obspass=observations(record,'PASS')
                for policy,olist in [('ALL',obsall),('PASS',obspass)]:
                    for obs in olist:
                        for scope in scopes:counts[scope,policy,obs['vtype'],obs['side']+'.'+obs['BD']]+=1
                passed=passes_filter(record)
                for obs in obsall:
                    if obs['vtype']!=a.vtype:continue
                    oid=f'r{n}_{obs["side"]}'
                    alleles,fallback=selected_alleles(record,obs['side'])
                    raw=[]
                    for ai,alt in alleles:
                        allele_id=f'{oid}a{ai}'
                        assert allele_id not in input_allele_ids
                        input_allele_ids.add(allele_id)
                        raw.append({'allele_index':ai,'start0':record.start,'ref':record.ref,'alt':alt})
                        af.write(f'{record.chrom}\t{record.pos}\t{oid}a{ai}\t{record.ref}\t{alt}\t.\tPASS\t.\n')
                    if fallback:qc['no_carried_gt_'+obs['side']+'_'+obs['BD']]+=1
                    if not raw:qc['observations_without_supported_alleles']+=1
                    d=dict(obs,oid=oid,chrom=record.chrom,pos1=record.pos,ref=record.ref,alts=list(record.alts or ()),genotype=list(record.samples[obs['side']].get('GT') or ()),gt_allele_fallback=fallback,filter=';'.join(record.filter.keys()) or '.',passes_filter=passed,scopes=sorted(scopes),raw_alleles=raw,BS=record.info.get('BS'),BK=record.samples[obs['side']].get('BK'),BI=record.samples[obs['side']].get('BI'))
                    mf.write(json.dumps(d)+'\n');qc['official_observations_selected']+=1
    comparison=compare_counts(counts,str(qprefix)+'.extended.csv',SCOPES)
    if not comparison['matched']:
        Path(str(prefix)+'.validation_failed.json').write_text(json.dumps(comparison,indent=2))
        raise RuntimeError('Official benchmark count conservation failed before distance assignment')
    normalized=work/'alleles.normalized.vcf.gz'
    command=[NORMALIZER,'norm','-f',cfg['reference'],'-c','e','-Oz','-o',str(normalized),str(allelevcf)]
    print(shlex.join(command),flush=True)
    with open(work/'normalization.log','w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=log)
    allelemap=defaultdict(list);normalized_allele_ids=set()
    with pysam.VariantFile(str(normalized)) as source:
        for record in source:
            assert record.id in input_allele_ids and record.id not in normalized_allele_ids,record.id
            normalized_allele_ids.add(record.id)
            oid=record.id.rsplit('a',1)[0]
            allelemap[oid].append({'pos1':record.pos,'start0':record.start,'ref':record.ref,'alt':record.alts[0]})
    assert normalized_allele_ids==input_allele_ids, 'Normalizer lost or added allele IDs'
    qc['allele_ids_conserved']=len(input_allele_ids)
    psv=PSVIndex(Path(a.prepared)/'psv_index.sqlite',{'chr'+str(i) for i in range(1,23)})
    scopeidx={name:BedIndex(read_bed(Path(a.prepared)/f'{name}.bed')) for name in SCOPES}
    familycounts=defaultdict(Counter);nonoverlap_counts=defaultdict(Counter);unambiguous_counts=defaultdict(Counter);normalized_overlap_counts=defaultdict(Counter)
    grouped=Counter();exclusions=Counter();normalization_changes=Counter();status=Counter();manual=[];rng=random.Random(20260929);seen_manual=0
    columns=['callset','tech','caller','vtype','scope','observation_id','side','all_decision','pass_decision','count','chrom','vcf_pos1','ref','alts','GT','FILTER','BS','BK','BI','genotype_allele_fallback','normalized_alleles','effective_spans0','state','distance_bin','distance_bp','family','copy_id','copy_start0','copy_end0','nearest_locus_id','nearest_psv_hash','nearest_psv_pos1','nearest_psv_ref','nearest_psv_alt','nearest_psv_effective_start0','nearest_psv_effective_end0','nearest_psv_effective_spans0','nearest_psv_type','nearest_partner_chrom','nearest_partner_pos1','nearest_partner_strand','nearest_nonoverlap_distance_bp','nearest_nonoverlap_locus_id','boundary_distance_bp','boundary_uncertain','raw_distance_bp','raw_distance_bin','normalization_changed_distance_bin','span_fully_within_scope','scope_overlap_effective_span','insertion_anchor','mixed_allele_types']
    with gzip.open(str(prefix)+'.sites.tsv.gz','wt') as outfh,gzip.open(meta,'rt') as inf:
        writer=csv.DictWriter(outfh,fieldnames=columns,delimiter='\t',extrasaction='ignore');writer.writeheader()
        for line in inf:
            obs=json.loads(line);alleles=allelemap.pop(obs['oid'],[])
            allele_runs=[allele_effective_spans(x['start0'],x['ref'],x['alt']) for x in alleles]
            spans_with_type=[span for runs in allele_runs for span in runs]
            if any(len(x['ref'])==len(x['alt']) and sum(r!=v for r,v in zip(x['ref'],x['alt']))>1 for x in alleles):
                qc['observations_with_multibase_substitutions']+=1
            if any(len(runs)>1 for runs in allele_runs):qc['observations_with_discontiguous_substitutions']+=1
            spans=[x[:2] for x in spans_with_type];kinds={x[2] for x in spans_with_type}
            if len(kinds)>1:qc['observations_mixed_allele_types']+=1
            annotation=psv.annotate(obs['chrom'],spans)
            rawspans=[span[:2] for x in obs['raw_alleles'] for span in allele_effective_spans(x['start0'],x['ref'],x['alt'])]
            rawannotation=psv.annotate(obs['chrom'],rawspans) if rawspans!=spans else annotation
            changed=annotation.get('distance_bin')!=rawannotation.get('distance_bin')
            status[annotation['state']]+=1
            normalization_changes['observations']+=1;normalization_changes['distance_bin_changed']+=int(changed)
            family=annotation.get('family','__'+annotation['state']+'__')
            locus=psv.locus(annotation['nearest_locus_id']) if annotation.get('nearest_locus_id') is not None else {}
            passdecision=('FN' if obs['BD']=='TP' and not obs['passes_filter'] else obs['BD']) if obs['side']=='TRUTH' else (obs['BD'] if obs['passes_filter'] else 'EXCLUDED_FILTER')
            for scope in obs['scopes']:
                normalized_scope_overlap=any(scopeidx[scope].overlaps(obs['chrom'],s,e) for s,e in spans)
                for policy,decision in [('ALL',obs['BD']),('PASS',passdecision)]:
                    klass=obs['side']+'.'+decision if decision!='EXCLUDED_FILTER' else 'QUERY.FILTERED'
                    key=(scope,policy,annotation['distance_bin'],family)
                    familycounts[key][klass]+=1
                    if normalized_scope_overlap:normalized_overlap_counts[key][klass]+=1
                    if klass!='QUERY.FILTERED':grouped[scope,policy,a.vtype,klass]+=1
                    if decision in ('UNK','EXCLUDED_FILTER'):exclusions[scope,policy,annotation['distance_bin'],decision]+=1
                    if annotation['state']=='usable':
                        if annotation.get('nearest_nonoverlap_distance_bp') is not None:
                            from psv_distance import distance_bin
                            b=distance_bin(annotation['nearest_nonoverlap_distance_bp']);nonoverlap_counts[scope,policy,b,family][klass]+=1
                        if not annotation.get('boundary_uncertain'):
                            unambiguous_counts[key][klass]+=1
                base={k:callset[k] for k in ('tech','caller')}
                base.update(annotation,callset=ident,vtype=a.vtype,scope=scope,observation_id=obs['oid'],side=obs['side'],all_decision=obs['BD'],pass_decision=passdecision,count=1,chrom=obs['chrom'],vcf_pos1=obs['pos1'],ref=obs['ref'],alts=','.join(obs['alts']),GT='/'.join('.' if g is None else str(g) for g in obs['genotype']),FILTER=obs['filter'],BS=obs['BS'],BK=obs['BK'],BI=obs['BI'],genotype_allele_fallback=obs['gt_allele_fallback'],normalized_alleles=json.dumps(alleles,separators=(',',':')),effective_spans0=json.dumps(spans,separators=(',',':')),family=family,nearest_psv_hash=locus.get('psv_hash'),nearest_psv_pos1=locus.get('raw_start0',-1)+1 if locus else None,nearest_psv_ref=locus.get('ref'),nearest_psv_alt=locus.get('alt'),nearest_psv_effective_start0=locus.get('start0'),nearest_psv_effective_end0=locus.get('end0'),nearest_psv_effective_spans0=json.dumps([x[:2] for x in allele_effective_spans(locus['raw_start0'],locus['ref'],locus['alt'])],separators=(',',':')) if locus else None,nearest_psv_type=locus.get('variant_class'),nearest_partner_chrom=locus.get('partner_chrom'),nearest_partner_pos1=locus.get('partner_start0',-1)+1 if locus else None,nearest_partner_strand=locus.get('partner_strand'),raw_distance_bp=rawannotation.get('distance_bp'),raw_distance_bin=rawannotation.get('distance_bin'),normalization_changed_distance_bin=changed,span_fully_within_scope=all(scopeidx[scope].contains(obs['chrom'],s,e) for s,e in spans) if spans else False,scope_overlap_effective_span=normalized_scope_overlap,insertion_anchor='INS_ANCHOR' in kinds,mixed_allele_types=len(kinds)>1)
                writer.writerow(base)
                qc['scope_rows_written']+=1
                if not base['span_fully_within_scope']:qc['scope_rows_not_fully_inside_literal_bed']+=1
                if not base['scope_overlap_effective_span']:qc['scope_rows_no_normalized_side_span_overlap']+=1
                if annotation.get('boundary_uncertain'):qc['scope_rows_boundary_uncertain']+=1
            if annotation['state']=='usable':
                seen_manual+=1
                example={'observation_id':obs['oid'],'side':obs['side'],'vtype':a.vtype,'chrom':obs['chrom'],'raw_pos1':obs['pos1'],'spans0':spans,**annotation,'nearest_psv':locus}
                if len(manual)<a.manual_checks:manual.append(example)
                else:
                    k=rng.randrange(seen_manual)
                    if k<a.manual_checks:manual[k]=example
    assert not allelemap, 'Unconsumed normalized observations'
    for key,value in counts.items():
        if key[2]==a.vtype:assert grouped[key]==value,(key,grouped[key],value)
    for ex in manual:
        rows=psv.psvis(ex['copy_id']).rows
        brute=min(span_distance(span,r[:2]) for span in ex['spans0'] for r in rows)
        assert brute==ex['distance_bp'],ex
        assert psv.copies[ex['copy_id']]['chrom']==ex['nearest_psv']['chrom']
        assert ex['nearest_psv']['family']==ex['family']
        ex['independent_bruteforce_distance']=brute
    write_counts(str(prefix)+'.family_counts.tsv',familycounts,callset,a.vtype)
    write_counts(str(prefix)+'.nonoverlap_psv_counts.tsv',nonoverlap_counts,callset,a.vtype)
    write_counts(str(prefix)+'.boundary_certain_counts.tsv',unambiguous_counts,callset,a.vtype)
    write_counts(str(prefix)+'.normalized_span_overlap_counts.tsv',normalized_overlap_counts,callset,a.vtype)
    with open(str(prefix)+'.exclusions.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['scope','policy','distance_bin','reason','n'])
        for k,n in sorted(exclusions.items()):w.writerow([*k,n])
    Path(str(prefix)+'.manual_checks.json').write_text(json.dumps(manual,indent=2)+'\n')
    report={'callset':ident,'vtype':a.vtype,'official_quantification_parity':comparison,'grouped_counts_conserved':True,'side_observation_unit':'official hap.py BVT/BD observation, never raw caller VCF row','qc':dict(qc),'state_counts':dict(status),'normalization_sensitivity':dict(normalization_changes),'distance_checks':len(manual),'normalization_command':command,'source_vcf':str(vcf.resolve()),'source_extended':str(Path(str(qprefix)+'.extended.csv').resolve()),'psv_index':str((Path(a.prepared)/'psv_index.sqlite').resolve()),'elapsed_seconds':time.time()-start,'job_id':os.environ.get('SLURM_JOB_ID'),'python':__import__('sys').version,'pysam':pysam.__version__}
    Path(str(prefix)+'.validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
