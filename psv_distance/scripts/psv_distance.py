"""Family- and local-copy-aware reference PSV interval distance queries."""
from bisect import bisect_left,bisect_right
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
import sqlite3
from prepare_psv_database import effective_span
from intervals import merge

DISTANCE_BINS=['0','1-50','51-150','151-300','301-1000','1001-5000','>5000']
SPECIAL_BINS=['no_usable_psv','ambiguous_copy','database_uncovered','representation_crosses_copy','unsupported_allele']

def distance_bin(d):
    if d is None:return 'no_usable_psv'
    for maximum,label in zip([0,50,150,300,1000,5000],DISTANCE_BINS):
        if d<=maximum:return label
    return '>5000'


def allele_effective_spans(start0,ref,alt):
    """Changed bases for substitutions; reference span/anchor for length changes.

    Equal-length alleles can have unchanged bases between substitutions. They
    contribute separate geometric runs but remain ONE official side observation.
    Identical REF/ALT yields no changed span and is not an insertion.
    """
    if len(ref)!=len(alt):return [effective_span(start0,ref,alt)]
    spans=[];run=None
    for i,(r,a) in enumerate(zip(ref.upper(),alt.upper())):
        if r!=a and run is None:run=i
        if r==a and run is not None:
            spans.append((start0+run,start0+i,'SNP' if i-run==1 else 'MNV'));run=None
    if run is not None:
        spans.append((start0+run,start0+len(ref),'SNP' if len(ref)-run==1 else 'MNV'))
    return spans

def span_distance(a,b):
    """Distance between reference BASES, not BED gap length: adjacent bases = 1."""
    s,e=a;t,u=b
    if e<=s or u<=t:raise ValueError('Empty interval; encode insertion as documented anchor base')
    return max(t-(e-1),s-(u-1),0)

class IntervalPSVs:
    def __init__(self,rows):
        # rows(start0,end0,locus_id). Repeated partner edges keep IDs, not extra distance weight.
        self.rows=sorted(rows)
        self.starts=[r[0] for r in self.rows]
        seen=set();self.multi_run_ids=set()
        for row in self.rows:
            if row[2] in seen:self.multi_run_ids.add(row[2])
            seen.add(row[2])
        self.by_end=sorted(self.rows,key=lambda r:(r[1],r[0],r[2]))
        self.ends=[r[1] for r in self.by_end]
        self.prefix_best=[];best=None
        for row in self.rows:
            if best is None or row[1]>best[1]:best=row
            self.prefix_best.append(best)
    def nearest(self,s,e,exclude_overlaps=False,excluded_loci=None):
        if not self.rows:return None,None
        excluded_loci=excluded_loci or set()
        if excluded_loci and not exclude_overlaps:
            raise ValueError('Locus exclusions require nonoverlap mode')
        if not exclude_overlaps:
            i=bisect_left(self.starts,e)-1
            if i>=0 and self.prefix_best[i][1]>s:
                return 0,self.prefix_best[i][2]
        candidates=[]
        i=bisect_right(self.ends,s)-1
        while i>=0 and self.by_end[i][2] in excluded_loci:i-=1
        if i>=0:
            row=self.by_end[i];candidates.append((span_distance((s,e),row[:2]),row[2]))
        i=bisect_left(self.starts,e)
        while i<len(self.rows) and self.rows[i][2] in excluded_loci:i+=1
        if i<len(self.rows):
            row=self.rows[i];candidates.append((span_distance((s,e),row[:2]),row[2]))
        return min(candidates) if candidates else (None,None)


    def overlapping_loci(self,s,e):
        """Parent PSV IDs touching any changed base in the query span."""
        i=bisect_left(self.starts,e)-1;out=set()
        while i>=0 and self.prefix_best[i][1]>s:
            row=self.rows[i]
            if row[1]>s:out.add(row[2])
            i-=1
        return out

class PSVIndex:
    def __init__(self,path,chroms=None):
        self.connection=sqlite3.connect('file:'+str(Path(path).resolve())+'?mode=ro',uri=True)
        self.connection.row_factory=sqlite3.Row
        self.connection.execute('PRAGMA cache_size=-65536')
        self.copies={r['copy_id']:dict(r) for r in self.connection.execute('SELECT * FROM copies') if chroms is None or r['chrom'] in chroms}
        self.by_chrom=defaultdict(list)
        for r in self.copies.values():self.by_chrom[r['chrom']].append(r)
        self.starts={};self.maxends={}
        for c,rows in self.by_chrom.items():
            rows.sort(key=lambda r:r['start0']);self.starts[c]=[r['start0'] for r in rows]
            ends=[];mx=0
            for r in rows:mx=max(mx,r['end0']);ends.append(mx)
            self.maxends[c]=ends
        self._sets={}
    def matching_copies(self,chrom,s,e):
        rows=self.by_chrom.get(chrom,[])
        if not rows:return []
        i=bisect_left(self.starts[chrom],e)-1;out=[]
        while i>=0 and self.maxends[chrom][i]>s:
            r=rows[i]
            if r['start0']<=s and r['end0']>=e:out.append(r['copy_id'])
            i-=1
        return out
    def overlapping_copies(self,chrom,s,e):
        """Find all physical copy envelopes touched by a reference span."""
        rows=self.by_chrom.get(chrom,[])
        if not rows:return []
        i=bisect_left(self.starts[chrom],e)-1;out=[]
        while i>=0 and self.maxends[chrom][i]>s:
            r=rows[i]
            if r['start0']<e and r['end0']>s:out.append(r['copy_id'])
            i-=1
        return out
    def psvis(self,copy_id):
        if copy_id not in self._sets:
            parents=self.connection.execute('SELECT raw_start0,ref,alt,locus_id FROM psv_loci WHERE copy_id=? AND usable=1 ORDER BY start0',(copy_id,))
            rows=[]
            for parent in parents:
                if parent['raw_start0'] is None or not parent['ref'] or not parent['alt']:
                    raise ValueError('Usable PSV lacks source alleles: '+str(parent['locus_id']))
                rows.extend((s,e,parent['locus_id']) for s,e,_ in allele_effective_spans(parent['raw_start0'],parent['ref'],parent['alt']))
            self._sets[copy_id]=IntervalPSVs(rows)
        return self._sets[copy_id]
    @lru_cache(maxsize=100000)
    def locus(self,locus_id):
        return dict(self.connection.execute('SELECT * FROM psv_loci WHERE locus_id=?',(int(locus_id),)).fetchone())
    def annotate(self,chrom,spans):
        if not spans:return {'state':'unsupported_allele','distance_bin':'unsupported_allele'}
        hits=[set(self.matching_copies(chrom,s,e)) for s,e in spans]
        shared=set.intersection(*hits)
        touched=set().union(*(set(self.overlapping_copies(chrom,s,e)) for s,e in spans))
        if not shared:
            status='representation_crosses_copy' if touched else 'database_uncovered'
            return {'state':status,'distance_bin':status,'candidate_copy_ids':';'.join(sorted(touched))}
        # A span partly entering another copy/family is conservatively ambiguous,
        # even if one envelope contains the entire representation.
        if len(shared)>1 or touched-shared:
            return {'state':'ambiguous_copy','distance_bin':'ambiguous_copy',
                    'candidate_copy_ids':';'.join(sorted(touched))}
        copy_id=next(iter(shared));copy=self.copies[copy_id];psvs=self.psvis(copy_id)
        result={'copy_id':copy_id,'family':copy['family'],'copy_start0':copy['start0'],'copy_end0':copy['end0']}
        candidates=[psvs.nearest(s,e) for s,e in spans];candidates=[x for x in candidates if x[0] is not None]
        if not candidates:
            return dict(result,state='no_usable_psv',distance_bin='no_usable_psv',boundary_uncertain=True)
        d,psvid=min(candidates)
        # Nearest PSV not overlapping ANY genotype-carried allele span.
        excluded_spans=merge(spans)
        if d!=0:
            nonself=(d,psvid)
        elif len(excluded_spans)==1:
            # Standard single-run PSV needs only interval exclusion. Multi-run
            # PSV must be removed as a whole parent, including distant runs.
            blocked=psvs.overlapping_loci(*excluded_spans[0]) & psvs.multi_run_ids if psvs.multi_run_ids else set()
            nonself=psvs.nearest(*excluded_spans[0],exclude_overlaps=True,excluded_loci=blocked)
        else:
            blocked=set().union(*(psvs.overlapping_loci(s,e) for s,e in excluded_spans))
            choices=[psvs.nearest(s,e,exclude_overlaps=True,excluded_loci=blocked) for s,e in excluded_spans]
            choices=[x for x in choices if x[0] is not None]
            nonself=min(choices) if choices else (None,None)
        edge=min(min(s-copy['start0']+1,copy['end0']-e+1) for s,e in spans)
        return dict(result,state='usable',distance_bp=d,distance_bin=distance_bin(d),nearest_locus_id=psvid,
                    nearest_nonoverlap_distance_bp=nonself[0],nearest_nonoverlap_locus_id=nonself[1],
                    boundary_distance_bp=edge,boundary_uncertain=d>=edge)
