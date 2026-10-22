"""Algorithmically distinct LV support-hypergraph orbit count via decorated core forests.

Only the Python standard library is used. This deliberately does not enumerate
cone spanning-tree bases or incidence matrices. It implements the leaf-quotient
bijection of Theorem 6.3 (core-forest classification) in the paper. Trees are generated independently by adding
one leaf and canonically deduplicating after each order.
"""
import argparse
import itertools
import json
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent.parent

def root_code(T,v,parent=-1):
    return '('+''.join(sorted(root_code(T,w,v) for w in T[v] if w!=parent))+')'

def canonical(T):
    return min(root_code(T,v) for v in range(len(T)))

def trees_up_to(nmax):
    layer=[(frozenset(),)]
    yield 1,layer
    for n in range(2,nmax+1):
        new={}
        for T in layer:
            for v in range(n-1):
                S=[set(a) for a in T]+[{v}]
                S[v].add(n-1)
                S=tuple(frozenset(a) for a in S)
                new[canonical(S)]=S
        layer=[new[k] for k in sorted(new)]
        yield n,layer

def edges(T):
    return [(v,w) for v in range(len(T)) for w in T[v] if v<w]

def auts(T,labels):
    groups={}
    for v in range(len(T)):groups.setdefault((len(T[v]),labels[v]),[]).append(v)
    groups=list(groups.values())
    es=set(edges(T));ans=[]
    for ps in itertools.product(*(itertools.permutations(g) for g in groups)):
        p=list(range(len(T)))
        for g,h in zip(groups,ps):
            for u,v in zip(g,h):p[u]=v
        if {tuple(sorted((p[u],p[v]))) for u,v in es}==es:ans.append(p)
    return ans

def core_count(T):
    n=len(T)
    if n<=2:return 1,1
    core=[v for v in range(n) if len(T[v])>=2]
    pos={v:i for i,v in enumerate(core)}
    H=tuple(frozenset(pos[w] for w in T[v] if w in pos) for v in core)
    leaves=[sum(len(T[w])==1 for w in T[v]) for v in core]
    autos=auts(H,leaves)
    es=edges(H)
    objects=set();raw=0
    for bits in range(1<<len(es)):
        selected=[e for k,e in enumerate(es) if bits>>k&1]
        F=[set() for _ in core]
        for u,v in selected:F[u].add(v);F[v].add(u)
        unused=set(range(len(core)));comps=[]
        while unused:
            todo=[unused.pop()];C=set(todo)
            while todo:
                v=todo.pop()
                for w in F[v]&unused:unused.remove(w);C.add(w);todo.append(w)
            comps.append(C)
        choices=[[(v,t) for v in sorted(C) for t in range(1+(leaves[v]>0))] for C in comps]
        for roots in itertools.product(*choices):
            raw+=1;images=[]
            for p in autos:
                eimage=tuple(sorted(tuple(sorted((p[u],p[v]))) for u,v in selected))
                rimage=tuple(sorted((p[v],t) for v,t in roots))
                images.append((eimage,rimage))
            objects.add(min(images))
    return len(objects),raw

def fib(n):
    a,b=0,1
    for _ in range(n):a,b=b,a+b
    return a

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',default=str(PACKAGE_ROOT/'data'/'lv_core_counts.json'))
    args=ap.parse_args()
    all_rows=[]
    expected_trees={1:1,2:1,3:1,4:2,5:3,6:6,7:11,8:23,9:47,10:106}
    for n,Ts in trees_up_to(10):
        assert len(Ts)==expected_trees[n]
        if n==1:continue
        rows=[]
        for T in Ts:
            q,raw=core_count(T)
            rows.append({'edges':edges(T),'q':q,'core_objects':raw})
            if max(map(len,T))<=2:assert q==(fib(2*n-2)+fib(n-1))//2
            if max(map(len,T))==n-1:assert q==(1 if n==2 else 2)
        print(json.dumps({'n':n,'trees':len(rows),'counts_sorted':sorted(x['q'] for x in rows)}),flush=True)
        all_rows.append({'n':n,'label':'exact computation','status':'PASS','trees':rows})
    with Path(args.out).open('w') as f:
        json.dump({'status':'PASS','by_n':all_rows},f,indent=2,sort_keys=True)
        f.write('\n')

if __name__=='__main__':main()
