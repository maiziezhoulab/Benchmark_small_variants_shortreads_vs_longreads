#!/usr/bin/env python3
"""Construct only requested evaluation scopes and common reference-base distance bands."""
import argparse,csv,json,sqlite3
from collections import defaultdict
from pathlib import Path
from intervals import *
from psv_distance import DISTANCE_BINS, PSVIndex

def unique_coverage(rows):
    events=defaultdict(lambda:defaultdict(list))
    for chrom,s,e,family,copy in rows:
        events[chrom][s].append((copy,1));events[chrom][e].append((copy,-1))
    unique=defaultdict(list);ambiguous=defaultdict(list);bycopy=defaultdict(list)
    for chrom,ev in events.items():
        active=defaultdict(int);previous=None
        for pos,changes in sorted(ev.items()):
            owners=[c for c,n in active.items() if n>0]
            if previous is not None and previous<pos:
                if len(owners)==1:
                    unique[chrom].append((previous,pos));bycopy[owners[0]].append((previous,pos))
                elif len(owners)>1:ambiguous[chrom].append((previous,pos))
            for c,n in changes:
                active[c]+=n
                if active[c]==0:del active[c]
            previous=pos
    return {c:merge(v) for c,v in unique.items()},{c:merge(v) for c,v in ambiguous.items()},{c:merge(v) for c,v in bycopy.items()}

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',default='config/inputs.json');ap.add_argument('--prepared',default='prepared');a=ap.parse_args()
    root=Path(a.prepared);cfg=json.load(open(a.config));db=sqlite3.connect(str(root/'psv_index.sqlite'))
    auto={'chr'+str(i) for i in range(1,23)}
    hc={c:v for c,v in read_bed(cfg['highconf_bed']).items() if c in auto};sd=read_bed(cfg['segdup_bed']);lm=read_bed(cfg['lowmap_bed'])
    initial=genomic_intersect(hc,sd)
    rows=[r for r in db.execute('SELECT chrom,start0,end0,family,copy_id FROM coverage') if r[0] in auto]
    raw=defaultdict(list)
    for c,s,e,fam,cp in db.execute('SELECT chrom,start0,end0,family,copy_id FROM raw_coverage'):
        if c in auto:raw[c].append((s,e))
    raw={c:merge(v) for c,v in raw.items()}
    unique,ambig,bycopy=unique_coverage(rows)
    allcov=defaultdict(list)
    for c,s,e,_,_ in rows:allcov[c].append((s,e))
    allcov={c:merge(v) for c,v in allcov.items()}
    scopes={'Segdup':genomic_intersect(initial,unique)}
    scopes['Segdup_lowmap']=genomic_intersect(scopes['Segdup'],lm)
    exclusions={'database_uncovered':genomic_subtract(initial,raw),
                'incompatible_or_unclear_homology':genomic_intersect(initial,genomic_subtract(raw,allcov)),
                'ambiguous_copy':genomic_intersect(initial,ambig)}
    for label,b in scopes.items():write_bed(root/f'{label}.bed',b)
    for label,b in exclusions.items():write_bed(root/f'excluded_{label}.bed',b)
    with open(root/'scopes.tsv','w') as f:
        for label in scopes:f.write(f'{label}\t{(root/f"{label}.bed").resolve()}\n')
    copies={r[0]:r for r in db.execute('SELECT copy_id,family,chrom,start0,end0 FROM copies')}
    bands=defaultdict(lambda:defaultdict(list));copy_support=[]
    psvindex=PSVIndex(root/'psv_index.sqlite',auto)
    for copy,cov in bycopy.items():
        _,family,chrom,cs,ce=copies[copy]
        scoped=intersect(cov,scopes['Segdup'].get(chrom,[]))
        if not scoped:continue
        psvs=merge([r[:2] for r in psvindex.psvis(copy).rows])
        copy_support.append((copy,family,chrom,cs,ce,len(psvs),sum(e-s for s,e in scoped)))
        if not psvs:bands['no_usable_psv'][chrom].extend(scoped);continue
        previous=[]
        for radius,label in zip([0,50,150,300,1000,5000],DISTANCE_BINS):
            expanded=merge([(max(cs,s-radius),min(ce,e+radius)) for s,e in psvs])
            current=intersect(scoped,expanded)
            bands[label][chrom].extend(subtract(current,previous));previous=current
        bands['>5000'][chrom].extend(subtract(scoped,previous))
    area_rows=[]
    for scope,scopebed in scopes.items():
        for label in DISTANCE_BINS+['no_usable_psv']:
            b=genomic_intersect(bands[label],scopebed)
            area_rows.append({'scope':scope,'distance_bin':label,'evaluable_bp':total_length(b)})
            write_bed(root/f'{scope}_distance_{label.replace(">","gt")}.bed',b)
        assert sum(r['evaluable_bp'] for r in area_rows if r['scope']==scope)==total_length(scopebed)
    with open(root/'scope_distance_areas.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=['scope','distance_bin','evaluable_bp'],delimiter='\t');w.writeheader();w.writerows(area_rows)
    with open(root/'copy_support.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['copy_id','family','chrom','start0','end0','n_psv_union_spans','evaluable_bp']);w.writerows(copy_support)
    report={'baseline_giab_segdup_bp':total_length(initial),'scopes_bp':{k:total_length(v) for k,v in scopes.items()},'excluded_bp':{k:total_length(v) for k,v in exclusions.items()},'copy_envelopes_in_scope':len(copy_support),'families_in_scope':len({r[1] for r in copy_support}),'scope_definition':'GIAB literal confident BED intersect segdup intersect reference-compatible unambiguous database homology windows; second scope additionally lowmap','area_definition':'Reference base distances; a multibase indel assigned by its effective span can differ from its anchor-base area bin','shared_bins':DISTANCE_BINS,'psv_geometry':'Equal-length alleles use contiguous actual mismatch runs; length-changing alleles use trimmed reference spans; insertions use left anchors'}
    assert total_length(scopes['Segdup'])+sum(total_length(b) for b in exclusions.values())==total_length(initial)
    (root/'scope_qc.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
