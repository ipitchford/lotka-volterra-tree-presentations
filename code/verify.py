#!/usr/bin/env python3
"""Exact, standard-library checks for generic Lotka--Volterra tree rigidity.

The all-dimensional proofs are in the paper. This script independently checks
coefficient identities, the completeness of linear Darboux supports at exact
sample parameters, spanning-tree counts, and support-hypergraph orbit counts.
Run from the package root: python3 code/verify.py --max-n 9 --out data/verification.json
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations
import json
from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parent.parent


def adjacency(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return [sorted(x) for x in adj]


def rooted_code(adj, v, parent=-1):
    return '(' + ''.join(sorted(rooted_code(adj, w, v)
                                for w in adj[v] if w != parent)) + ')'


def centres(adj):
    live = set(range(len(adj)))
    while len(live) > 2:
        leaves = {v for v in live if sum(w in live for w in adj[v]) <= 1}
        live -= leaves
    return sorted(live)


def tree_code(n, edges):
    adj = adjacency(n, edges)
    return min(rooted_code(adj, r) for r in centres(adj))


def all_trees(max_n):
    result = {1: {'()': ()}}
    for n in range(2, max_n + 1):
        nxt = {}
        for edges in result[n - 1].values():
            for v in range(n - 1):
                grown = tuple(sorted((*edges, (v, n - 1))))
                nxt.setdefault(tree_code(n, grown), grown)
        result[n] = dict(sorted(nxt.items()))
    return result


def exact_det(matrix):
    """Bareiss elimination; integer input and exact integer divisions."""
    a = [list(map(int, row)) for row in matrix]
    n, sign, old = len(a), 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            r = next((r for r in range(k + 1, n) if a[r][k]), None)
            if r is None:
                return 0
            a[k], a[r] = a[r], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j] * pivot - a[i][k] * a[k][j]
                assert value % old == 0
                a[i][j] = value // old
            a[i][k] = 0
        old = pivot
    return sign * a[-1][-1] if n else 1


def matrix_tree_count(n, edges):
    a = [[int(i == j) for j in range(n)] for i in range(n)]
    for u, v in edges:
        a[u][u] += 1
        a[v][v] += 1
        a[u][v] -= 1
        a[v][u] -= 1
    return exact_det(a)


def rooted_iso(adj, v, pv, w, pw, mapping):
    mapping[v] = w
    children_v = sorted((rooted_code(adj, x, v), x) for x in adj[v] if x != pv)
    children_w = sorted((rooted_code(adj, x, w), x) for x in adj[w] if x != pw)
    assert [x[0] for x in children_v] == [x[0] for x in children_w]
    for (_, x), (_, y) in zip(children_v, children_w):
        rooted_iso(adj, x, v, y, w, mapping)


def tree_automorphism_generators(n, edges):
    adj = adjacency(n, edges)
    generators = []

    def visit(v, parent):
        groups = defaultdict(list)
        for w in adj[v]:
            if w != parent:
                groups[rooted_code(adj, w, v)].append(w)
                visit(w, v)
        for group in groups.values():
            for x, y in zip(group, group[1:]):
                mapping = {}
                rooted_iso(adj, x, v, y, v, mapping)
                p = list(range(n))
                for a, b in mapping.items():
                    p[a], p[b] = b, a
                generators.append(tuple(p))

    cs = centres(adj)
    visit(cs[0], -1)
    if len(cs) == 2:
        u, v = cs
        if rooted_code(adj, u, v) == rooted_code(adj, v, u):
            mapping = {}
            rooted_iso(adj, u, v, v, u, mapping)
            p = list(range(n))
            for a, b in mapping.items():
                p[a], p[b] = b, a
            generators.append(tuple(p))
    eset = set(edges)
    for p in generators:
        assert sorted(p) == list(range(n))
        assert {tuple(sorted((p[u], p[v]))) for u, v in edges} == eset
    return list(dict.fromkeys(generators))


def cone_edges(n, edges):
    # The cone vertex is n internally, and 0 in the manuscript.
    return [(v, n) for v in range(n)] + list(edges)


def matroid_generators(n, edges):
    gedges = cone_edges(n, edges)
    index = {tuple(sorted(e)): i for i, e in enumerate(gedges)}
    result = []
    if n == 2:
        return [(1, 0, 2), (0, 2, 1)]
    for p in tree_automorphism_generators(n, edges):
        p = (*p, n)
        result.append(tuple(index[tuple(sorted((p[u], p[v])))] for u, v in gedges))
    for v, nbrs in enumerate(adjacency(n, edges)):
        if len(nbrs) == 1:
            j = index[tuple(sorted((v, nbrs[0])))]
            p = list(range(2 * n - 1))
            p[v], p[j] = j, v
            result.append(tuple(p))
    return list(dict.fromkeys(result))


def bases(n, edges):
    gedges = cone_edges(n, edges)
    output = []
    for chosen in combinations(range(2 * n - 1), n):
        parent = list(range(n + 1))

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v

        for j in chosen:
            u, v = gedges[j]
            u, v = find(u), find(v)
            if u == v:
                break
            parent[u] = v
        else:
            output.append(sum(1 << j for j in chosen))
    return output


def apply_permutation(mask, p):
    result = 0
    while mask:
        bit = mask & -mask
        j = bit.bit_length() - 1
        result |= 1 << p[j]
        mask -= bit
    return result


def orbit_counts(n, edges, bs):
    generators = matroid_generators(n, edges)
    remaining = set(bs)
    sizes, representatives = [], []
    while remaining:
        seed = min(remaining)
        orbit = {seed}
        stack = [seed]
        while stack:
            current = stack.pop()
            for p in generators:
                image = apply_permutation(current, p)
                assert image in remaining or image in orbit
                if image not in orbit:
                    orbit.add(image)
                    stack.append(image)
        remaining -= orbit
        sizes.append(len(orbit))
        representatives.append(seed)
    assert sum(sizes) == len(bs)
    return sizes, representatives


def sample_matrix(n, edges, signed=False):
    adj = adjacency(n, edges)
    a = [j + 1 for j in range(n)]
    matrix = [[0] * n for _ in range(n)]
    for j in range(n):
        matrix[j][j] = a[j]
        for k, v in enumerate(adj[j]):
            offset = 2 ** (k + 1)
            if signed and k % 2:
                offset = -offset
            value = a[j] + offset
            stack = [(v, j)]
            while stack:
                x, px = stack.pop()
                matrix[x][j] = value
                stack.extend((y, x) for y in adj[x] if y != px)
    return matrix


def direct_darboux_supports(a):
    """Independent coefficient comparison over every nonempty support.

    Off-support cofactor equality is checked directly. Nonzero pair
    differences determine all coefficients from the first coefficient.
    No connected-support or matroid theorem is used in this test.
    """
    n = len(a)
    supports = []
    for size in range(1, n + 1):
        for I in combinations(range(n), size):
            i = I[0]
            outside = set(range(n)) - set(I)
            if any(a[j][k] != a[i][k] for j in I for k in outside):
                continue
            c = {i: Q(1)}
            for j in I[1:]:
                assert a[j][i] != a[i][i]
                c[j] = -Q(a[i][j] - a[j][j], a[j][i] - a[i][i])
            if any(c[j] == 0 for j in I):
                continue
            if any(c[j] * (a[j][k] - a[k][k]) + c[k] * (a[k][j] - a[j][j])
                   for j, k in combinations(I, 2)):
                continue
            supports.append(I)
            # Verify the complete quadratic identity, including every zero
            # coefficient. This calculation is separate from support screening.
            full_c = [c.get(j, Q(0)) for j in range(n)]
            cofactor = [a[j][j] if j in I else a[i][j] for j in range(n)]
            for j in range(n):
                assert full_c[j] * (a[j][j] - cofactor[j]) == 0
                for k in range(j + 1, n):
                    assert full_c[j] * a[j][k] + full_c[k] * a[k][j] == \
                        full_c[j] * cofactor[k] + full_c[k] * cofactor[j]
    return set(supports)


SOURCE_TABLE_7 = {
    2: [1], 3: [2], 4: [2, 5], 5: [2, 8, 12],
    6: [2, 8, 5, 18, 21, 30],
    7: [2, 8, 8, 18, 21, 30, 22, 79, 12, 55, 76],
    8: [2, 8, 8, 18, 21, 5, 30, 30, 32, 79, 21, 55, 18, 48, 60,
        79, 126, 79, 112, 207, 30, 144, 195],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--max-n', type=int, default=9)
    ap.add_argument('--out', default=str(PACKAGE_ROOT / 'data' / 'verification.json'))
    args = ap.parse_args()
    records, checks = [], []
    trees = all_trees(args.max_n)
    known_tree_counts = {2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47, 10: 106}
    for n in range(2, args.max_n + 1):
        if n in known_tree_counts:
            assert len(trees[n]) == known_tree_counts[n]
        orbit_multiset = []
        for code, edges in trees[n].items():
            expected = {(i,) for i in range(n)} | set(edges)
            for signed in (False, True):
                assert direct_darboux_supports(sample_matrix(n, edges, signed)) == expected
            bs = bases(n, edges)
            determinant = matrix_tree_count(n, edges)
            assert len(bs) == determinant
            sizes, reps = orbit_counts(n, edges, bs)
            orbit_multiset.append(len(sizes))
            records.append({'n': n, 'tree_code': code, 'edges_zero_based': edges,
                            'basis_count': len(bs), 'det_I_plus_L': determinant,
                            'hypergraph_orbit_count': len(sizes),
                            'orbit_sizes': sizes, 'orbit_representative_masks': reps})
        if n in SOURCE_TABLE_7:
            assert sorted(orbit_multiset) == sorted(SOURCE_TABLE_7[n]), (n, sorted(orbit_multiset))
            checks.append({'n': n, 'label': 'exact computation', 'status': 'PASS', 'source_table_7_multiset_matches': True,
                           'tree_count': len(trees[n]),
                           'hypergraph_total': sum(orbit_multiset),
                           'sorted_orbit_counts': sorted(orbit_multiset)})
        else:
            checks.append({'n': n, 'label': 'exact computation', 'status': 'PASS', 'tree_count': len(trees[n]),
                           'hypergraph_total': sum(orbit_multiset),
                           'sorted_orbit_counts': sorted(orbit_multiset)})
        print(f'n={n}: {len(trees[n])} trees; {sum(orbit_multiset)} hypergraph classes; all checks pass', flush=True)
    output = {'status': 'PASS', 'arithmetic': 'Python integers and Fraction; no floating-point proof checks',
              'max_n': args.max_n, 'tree_count': len(records),
              'parameter_samples_per_tree': 2,
              'support_subsets_examined': sum(2 * (2 ** r['n'] - 1) for r in records),
              'basis_count_total': sum(r['basis_count'] for r in records),
              'source_table': {'url': 'https://arxiv.org/html/2411.18264v4', 'locator': 'Table 7'},
              'summary': checks, 'trees': records}
    Path(args.out).write_text(json.dumps(output, indent=2, sort_keys=True) + '\n')
    print(f'PASS: {len(records)} trees; wrote {Path(args.out).name}', flush=True)


if __name__ == '__main__':
    main()
