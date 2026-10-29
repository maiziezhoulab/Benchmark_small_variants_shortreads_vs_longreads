"""Half-open interval operations; no BED row is treated as a biological copy."""
from bisect import bisect_left, bisect_right
from collections import defaultdict
import gzip


def merge(intervals):
    out=[]
    for s,e in sorted(intervals):
        if e<=s: continue
        if out and s<=out[-1][1]:out[-1]=(out[-1][0],max(e,out[-1][1]))
        else:out.append((s,e))
    return out


def intersect(a,b):
    a=merge(a);b=merge(b);out=[];i=j=0
    while i<len(a) and j<len(b):
        s=max(a[i][0],b[j][0]);e=min(a[i][1],b[j][1])
        if s<e:out.append((s,e))
        if a[i][1]<b[j][1]:i+=1
        else:j+=1
    return out


def subtract(a,b):
    a=merge(a);b=merge(b);out=[];j=0
    for s,e in a:
        while j<len(b) and b[j][1]<=s:j+=1
        k=j;p=s
        while k<len(b) and b[k][0]<e:
            if p<b[k][0]:out.append((p,min(e,b[k][0])))
            p=max(p,b[k][1]);k+=1
            if p>=e:break
        if p<e:out.append((p,e))
    return out


def read_bed(path):
    d=defaultdict(list)
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            x=line.split();d[x[0]].append((int(x[1]),int(x[2])))
    return {c:merge(v) for c,v in d.items()}


def write_bed(path,d):
    with open(path,'w') as f:
        for c in sorted(d):
            for s,e in merge(d[c]):f.write(f'{c}\t{s}\t{e}\n')


def genomic_intersect(a,b):
    return {c:intersect(v,b.get(c,[])) for c,v in a.items()}


def genomic_subtract(a,b):
    return {c:subtract(v,b.get(c,[])) for c,v in a.items()}


def total_length(d):return sum(e-s for vals in d.values() for s,e in vals)


class BedIndex:
    def __init__(self,d):
        self.data={c:merge(v) for c,v in d.items()}
        self.starts={c:[s for s,e in v] for c,v in self.data.items()}
    def overlaps(self,c,s,e):
        vals=self.data.get(c,[])
        if not vals:return False
        i=bisect_left(self.starts[c],e)-1
        return i>=0 and vals[i][1]>s
    def contains(self,c,s,e):
        vals=self.data.get(c,[])
        if not vals:return False
        i=bisect_right(self.starts[c],s)-1
        return i>=0 and vals[i][1]>=e
