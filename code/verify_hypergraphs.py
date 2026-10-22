"""Direct small hypergraph-isomorphism check, independent of matroid symmetries.

For every unlabelled tree of orders 2 through 6, compute the fundamental paths
for every cone spanning tree and canonicalise their supports under all n!
coordinate permutations. Compare with generator-orbit counts from verify.py.
No automorphism-group assertion is used to canonicalise a hypergraph.
Only the Python standard library is needed.
"""
import argparse
from itertools import permutations
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
PACKAGE_ROOT = Path(__file__).resolve().parent.parent
from verify import all_trees, bases, cone_edges, orbit_counts


def fundamental_hypergraph(n, edges, basis):
    graph_edges = cone_edges(n, edges)
    chosen = [i for i in range(len(graph_edges)) if (basis >> i) & 1]
    adjacency = [[] for _ in range(n + 1)]
    for coordinate, edge_index in enumerate(chosen):
        u, v = graph_edges[edge_index]
        adjacency[u].append((v, coordinate))
        adjacency[v].append((u, coordinate))
    supports = []
    for edge_index, (u, v) in enumerate(graph_edges):
        if (basis >> edge_index) & 1:
            continue
        stack = [(u, -1, 0)]
        while stack:
            w, parent, support = stack.pop()
            if w == v:
                supports.append(support)
                break
            stack.extend((t, w, support | (1 << c))
                         for t, c in adjacency[w] if t != parent)
        else:
            raise AssertionError('A basis did not connect an omitted edge.')
    assert len(supports) == len(set(supports)) == n - 1
    return tuple(sorted(supports))


def permutation_tables(n):
    tables = []
    for p in permutations(range(n)):
        tables.append(tuple(sum(1 << p[j] for j in range(n) if (mask >> j) & 1)
                            for mask in range(1 << n)))
    return tables


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=str(PACKAGE_ROOT / 'data' / 'direct_hypergraph_verification.json'))
    args = ap.parse_args()
    trees = all_trees(6)
    records = []
    for n in range(2, 7):
        tables = permutation_tables(n)
        for edges in trees[n].values():
            bs = list(bases(n, edges))
            raw_supports = {fundamental_hypergraph(n, edges, b) for b in bs}
            canonical = {min(tuple(sorted(table[s] for s in h)) for table in tables)
                         for h in raw_supports}
            orbit_sizes, _ = orbit_counts(n, edges, bs)
            assert len(canonical) == len(orbit_sizes)
            records.append({'n': n, 'label': 'exact computation', 'status': 'PASS', 'edges': edges, 'bases': len(bs),
                            'direct_hypergraph_classes': len(canonical),
                            'generator_orbit_classes': len(orbit_sizes)})
    result = {'status': 'PASS', 'max_n': 6, 'trees_checked': len(records),
              'bases_checked': sum(r['bases'] for r in records),
              'method': 'Fundamental paths, followed by all n! hypergraph vertex relabellings.',
              'records': records}
    Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'records'}))


if __name__ == '__main__':
    main()
