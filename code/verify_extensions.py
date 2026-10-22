#!/usr/bin/env python3
"""Exact checks for the extensions of the Lotka--Volterra tree classification.

Notation. T is a tree on V = {0, ..., n-1}. Parameters a_j and b_{v->j}, one
for every oriented edge v->j. The tree matrix has A_jj = a_j and
A_ij = b_{v(i,j)->j}, where v(i,j) is the neighbour of j on the path j..i.
d_{v->j} = b_{v->j} - a_j. Delta_T = prod d_{v->j} *
prod_j prod_{u<v in N(j)} (d_{u->j} - d_{v->j}) (d_{u->j} + d_{v->j}).
The edge ray of uv is P_uv = d_{v->u} x_u - d_{u->v} x_v.

All arithmetic is exact (Python integers and fractions.Fraction); sympy is used
only for the symbolic checks (g) and (i). Pseudo-random points use fixed
string seeds, so the output is byte-identical between runs.

Run from the package root:  python3 code/verify_extensions.py
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
PACKAGE_ROOT = Path(__file__).resolve().parent.parent

import verify as V  # noqa: E402
import lv_core_count as C  # noqa: E402
from verify_hypergraphs import fundamental_hypergraph  # noqa: E402

EXACT = 'exact computation'
EXPLORATORY = 'exploratory'


# --------------------------------------------------------------------------
# Exact linear algebra over Q
# --------------------------------------------------------------------------

def rref(rows, ncols):
    m = [[Q(x) for x in r] for r in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        inv = 1 / m[r][c]
        m[r] = [x * inv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    return m[:r], pivots


def rank(rows, ncols=None):
    if not rows:
        return 0
    if ncols is None:
        ncols = len(rows[0])
    return len(rref(rows, ncols)[1])


def nullspace(rows, ncols):
    if not rows:
        return [[Q(int(i == j)) for j in range(ncols)] for i in range(ncols)]
    m, pivots = rref(rows, ncols)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for f in free:
        v = [Q(0)] * ncols
        v[f] = Q(1)
        for row, p in zip(m, pivots):
            v[p] = -row[f]
        basis.append(v)
    return basis


def inverse(mat):
    n = len(mat)
    aug = [list(map(Q, row)) + [Q(int(i == j)) for j in range(n)]
           for i, row in enumerate(mat)]
    m, pivots = rref(aug, 2 * n)
    assert pivots[:n] == list(range(n)), 'singular matrix'
    return [row[n:] for row in m]


def normalise(vec):
    first = next(x for x in vec if x != 0)
    return tuple(Q(x) / first for x in vec)


def fstr(x):
    x = Q(x)
    return str(x.numerator) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'


def vstr(vec):
    return [fstr(x) for x in vec]


# --------------------------------------------------------------------------
# Complete linear Darboux polynomials of a numeric LV matrix
# --------------------------------------------------------------------------

def complete_rays(A, r=None):
    """Every nonempty support I of x' = x * (r + A x), r constant vector.

    Coefficient comparison in P' = P * (C0 + C(x)), P = sum_{i in I} c_i x_i,
    all c_i nonzero, gives: r_i = C0 constant on I (linear terms); C_i = A_ii
    for i in I; A_ij = A_kj = C_j for i, k in I, j outside I; and
    c_i (A_ij - A_jj) + c_j (A_ji - A_ii) = 0 for i < j in I.
    Returns supports admitting a full-support solution with the solution
    dimension; dimension 1 supports give rays, dimension >= 2 pencils.
    """
    n = len(A)
    out = []
    for size in range(1, n + 1):
        for I in combinations(range(n), size):
            if r is not None and any(r[i] != r[I[0]] for i in I):
                continue
            Iset = set(I)
            if any(A[i][j] != A[I[0]][j] for i in I for j in range(n) if j not in Iset):
                continue
            pos = {v: k for k, v in enumerate(I)}
            rows = []
            for i, j in combinations(I, 2):
                row = [Q(0)] * size
                row[pos[i]] += A[i][j] - A[j][j]
                row[pos[j]] += A[j][i] - A[i][i]
                rows.append(row)
            basis = nullspace(rows, size)
            if not basis or any(all(b[k] == 0 for b in basis) for k in range(size)):
                continue
            rec = {'support': I, 'dim': len(basis)}
            if len(basis) == 1:
                full = [Q(0)] * n
                for v, k in pos.items():
                    full[v] = basis[0][k]
                rec['vector'] = normalise(full)
            rec['cofactor'] = tuple(Q(A[j][j]) if j in Iset else Q(A[I[0]][j])
                                    for j in range(n))
            out.append(rec)
    return out


def ray_set(cr):
    return {(r['support'], r['vector']) for r in cr if r['dim'] == 1}


def pencils(cr):
    return [r for r in cr if r['dim'] >= 2]


def darboux_identity_holds(A, c, cof, r=None, C0=0):
    """Direct check that P(x * (r + A x)) == P(x) (C0 + C(x)) coefficientwise."""
    n = len(A)
    for i in range(n):
        lin = c[i] * ((r[i] if r else 0) - C0)
        if lin != 0:
            return False
        if c[i] * (A[i][i] - cof[i]) != 0:
            return False
        for j in range(i + 1, n):
            if c[i] * A[i][j] + c[j] * A[j][i] != c[i] * cof[j] + c[j] * cof[i]:
                return False
    return True


# --------------------------------------------------------------------------
# Trees, tree matrices, rays of D_T
# --------------------------------------------------------------------------

def adjacency(n, edges):
    return V.adjacency(n, edges)


def first_step(adj, j):
    """first[i] = neighbour of j on the path j..i (None for i == j)."""
    first = {j: None}
    stack = [(w, w) for w in adj[j]]
    seen = {j}
    while stack:
        x, f = stack.pop()
        if x in seen:
            continue
        seen.add(x)
        first[x] = f
        stack.extend((y, f) for y in adj[x] if y not in seen)
    return first


def toward(adj, targets):
    """w[j] = neighbour of j towards the vertex set targets (j not in targets)."""
    w = {}
    frontier = sorted(targets)
    seen = set(targets)
    while frontier:
        nxt = []
        for x in frontier:
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    w[y] = x
                    nxt.append(y)
        frontier = nxt
    return w


def tree_matrix(n, edges, a, d):
    adj = adjacency(n, edges)
    A = [[Q(0)] * n for _ in range(n)]
    for j in range(n):
        first = first_step(adj, j)
        for i in range(n):
            A[i][j] = Q(a[j]) if i == j else Q(a[j] + d[(first[i], j)])
    return A


def oriented(edges):
    return sorted([(u, v) for u, v in edges] + [(v, u) for u, v in edges])


def generic_point(n, edges, rng, signed=False):
    adj = adjacency(n, edges)
    a = [rng.randint(-9, 9) for _ in range(n)]
    d = {}
    for j in range(n):
        vals = rng.sample(range(1, 60), len(adj[j]))
        for v, x in zip(adj[j], vals):
            d[(v, j)] = -x if signed and rng.random() < 0.5 else x
    return a, d


def delta_factors(n, edges, d):
    adj = adjacency(n, edges)
    out = [(('d', v, j), d[(v, j)]) for v, j in oriented(edges)]
    for j in range(n):
        for u, v in combinations(adj[j], 2):
            out.append((('minus', j, u, v), d[(u, j)] - d[(v, j)]))
            out.append((('plus', j, u, v), d[(u, j)] + d[(v, j)]))
    return out


def delta_value(n, edges, d):
    p = 1
    for _, val in delta_factors(n, edges, d):
        p *= val
    return p


def tree_rays(n, edges, d):
    """Labelled rays of D_T: spokes ('s', i) and tree edges ('r', u, v)."""
    rays = {}
    for i in range(n):
        vec = [Q(0)] * n
        vec[i] = Q(1)
        rays[('s', i)] = vec
    for u, v in edges:
        vec = [Q(0)] * n
        vec[u] = Q(d[(v, u)])
        vec[v] = Q(-d[(u, v)])
        rays[('r', u, v)] = vec
    return rays


def ray_vertices(label):
    return {label[1]} if label[0] == 's' else {label[1], label[2]}


def lemma_cofactor(n, edges, a, d, label):
    adj = adjacency(n, edges)
    Ve = ray_vertices(label)
    w = toward(adj, Ve)
    return tuple(Q(a[j]) if j in Ve else Q(a[j] + d[(w[j], j)]) for j in range(n))


def in_span(target, gens):
    return rank(gens + [target]) == rank(gens)


def expected_tree_rayset(n, edges, d):
    return {(tuple(i for i in range(n) if vec[i] != 0), normalise(vec))
            for vec in tree_rays(n, edges, d).values()}


def is_tree_form(n, edges, cr):
    if pencils(cr):
        return False
    supports = sorted(r['support'] for r in cr)
    want = sorted([(i,) for i in range(n)] + [tuple(sorted(e)) for e in edges])
    return supports == want


def label_str(label):
    return f's{label[1]}' if label[0] == 's' else f'r{label[1]}{label[2]}'


def trees_upto(max_n):
    return V.all_trees(max_n)


def seeded(*parts):
    return random.Random('|'.join(map(str, parts)))


def result(label, failures, **details):
    return {'label': label, 'status': 'PASS' if not failures else 'FAIL',
            'failures': failures[:20], 'failure_count': len(failures), **details}


# --------------------------------------------------------------------------
# (a) Cofactor lemma and the triangle criterion
# --------------------------------------------------------------------------

def compat_relation(n, edges, a, d, A):
    rays = tree_rays(n, edges, d)
    labels = sorted(rays)
    cof = {lab: lemma_cofactor(n, edges, a, d, lab) for lab in labels}
    compat = {}
    for e, f in combinations(labels, 2):
        diff = [x - y for x, y in zip(cof[e], cof[f])]
        compat[(e, f)] = in_span(diff, [rays[e], rays[f]])
    return rays, labels, cof, compat


def triangle_pair(e, f, edges):
    es = {tuple(sorted(x)) for x in edges}
    for g, h in ((e, f), (f, e)):
        if g[0] == 's' and h[0] == 's' and tuple(sorted((g[1], h[1]))) in es:
            return True
        if g[0] == 's' and h[0] == 'r' and g[1] in (h[1], h[2]):
            return True
    return False


def check_a():
    fails, per_n = [], {}
    pairs = rays_checked = 0
    for n in range(2, 8):
        cnt = {'trees': 0, 'points': 0, 'ray_pairs': 0, 'compatible_pairs': 0}
        for code, edges in trees_upto(n)[n].items():
            cnt['trees'] += 1
            for k in range(2):
                a, d = generic_point(n, edges, seeded('a', n, code, k))
                A = tree_matrix(n, edges, a, d)
                cnt['points'] += 1
                cr = complete_rays(A)
                if ray_set(cr) != expected_tree_rayset(n, edges, d) or pencils(cr):
                    fails.append(f'n={n} {code} point {k}: ray set is not D_T')
                rays, labels, cof, compat = compat_relation(n, edges, a, d, A)
                by_support = {r['support']: r for r in cr}
                for lab in labels:
                    rays_checked += 1
                    sup = tuple(sorted(ray_vertices(lab)))
                    rec = by_support.get(sup)
                    if rec is None or rec['cofactor'] != cof[lab]:
                        fails.append(f'n={n} {code}: cofactor lemma {label_str(lab)}')
                    if not darboux_identity_holds(A, rays[lab], cof[lab]):
                        fails.append(f'n={n} {code}: identity {label_str(lab)}')
                for (e, f), ok in compat.items():
                    pairs += 1
                    cnt['ray_pairs'] += 1
                    cnt['compatible_pairs'] += ok
                    if ok != triangle_pair(e, f, edges):
                        fails.append(f'n={n} {code}: pair {label_str(e)},{label_str(f)}')
        per_n[str(n)] = cnt
    return result(EXACT, fails, statement=(
        'At generic exact points the complete ray set is D_T; every ray cofactor '
        'equals the lemma formula and satisfies the Darboux identity; '
        'c(e)-c(f) lies in span(l_e,l_f) iff e,f lie in a common triangle {s_u,s_v,r_uv}.'),
        rays_checked=rays_checked, ray_pairs_checked=pairs, per_n=per_n)


# --------------------------------------------------------------------------
# (b) One-sided rigidity
# --------------------------------------------------------------------------

def prufer_trees(n):
    if n == 1:
        yield ()
        return
    if n == 2:
        yield ((0, 1),)
        return
    for seq in product(range(n), repeat=n - 2):
        degree = [1] * n
        for x in seq:
            degree[x] += 1
        es = []
        for x in seq:
            leaf = min(i for i in range(n) if degree[i] == 1)
            es.append(tuple(sorted((leaf, x))))
            degree[leaf] -= 1
            degree[x] -= 1
        u, w = [i for i in range(n) if degree[i] == 1]
        es.append((u, w))
        yield tuple(sorted(es))


def lt_equations(n, edges):
    """Linear equations on the n^2 entries (index i*n+j) defining L_T."""
    adj = adjacency(n, edges)
    rows = []
    for j in range(n):
        first = first_step(adj, j)
        for v in adj[j]:
            branch = sorted(i for i in range(n) if i != j and first[i] == v)
            for i in branch[1:]:
                row = [0] * (n * n)
                row[branch[0] * n + j] = 1
                row[i * n + j] = -1
                rows.append(row)
    return rows


def in_family(M, n, edges):
    flat = [M[i][j] for i in range(n) for j in range(n)]
    return all(sum(r * x for r, x in zip(row, flat)) == 0 for row in lt_equations(n, edges))


def check_b():
    fails = []
    per_n = {}
    for n in range(2, 7):
        cnt = {'trees': 0, 'points': 0, 'bases': 0, 'compat_graph_spanning_tree': 0,
               'compat_graph_isomorphic_to_T': 0, 'labelled_target_trees_tested': 0,
               'family_memberships': 0}
        labelled = list(prufer_trees(n)) if n <= 5 else None
        for code, edges in trees_upto(n)[n].items():
            cnt['trees'] += 1
            cone = V.cone_edges(n, edges)
            # ray label in cone-edge order: spoke (v, n) -> ('s', v); edge -> ('r', u, v)
            cone_labels = [('s', v) for v in range(n)] + [('r', u, v) for u, v in edges]
            basis_masks = set(V.bases(n, edges))
            for k in range(2):
                a, d = generic_point(n, edges, seeded('b', n, code, k))
                A = tree_matrix(n, edges, a, d)
                cnt['points'] += 1
                rays, labels, cof, compat = compat_relation(n, edges, a, d, A)
                vecs = [rays[lab] for lab in cone_labels]
                indep = {sum(1 << j for j in S) for S in combinations(range(2 * n - 1), n)
                         if rank([vecs[j] for j in S]) == n}
                if indep != basis_masks:
                    fails.append(f'n={n} {code}: ray matroid bases differ from cone spanning trees')
                for mask in sorted(basis_masks):
                    cnt['bases'] += 1
                    S = [j for j in range(2 * n - 1) if mask >> j & 1]
                    B = [vecs[j] for j in S]
                    Binv = inverse(B)
                    # fundamental supports, linearly, in z-coordinates
                    lin_sup = []
                    for g in range(2 * n - 1):
                        if mask >> g & 1:
                            continue
                        coeffs = [sum(vecs[g][t] * Binv[t][s] for t in range(n)) for s in range(n)]
                        lin_sup.append(sum(1 << s for s in range(n) if coeffs[s] != 0))
                    if tuple(sorted(lin_sup)) != fundamental_hypergraph(n, edges, mask):
                        fails.append(f'n={n} {code}: linear supports differ from fundamental paths')
                    size2 = {tuple(s for s in range(n) if m >> s & 1)
                             for m in lin_sup if bin(m).count('1') == 2}
                    cg = set()
                    for s, t in combinations(range(n), 2):
                        e, f = sorted((cone_labels[S[s]], cone_labels[S[t]]))
                        if compat[(e, f)]:
                            cg.add((s, t))
                    if cg != size2:
                        fails.append(f'n={n} {code}: compatibility graph != size-2 fundamental paths')
                    spanning = len(cg) == n - 1 and len(_components(n, cg)) == 1
                    if spanning:
                        cnt['compat_graph_spanning_tree'] += 1
                        if V.tree_code(n, tuple(sorted(cg))) == code:
                            cnt['compat_graph_isomorphic_to_T'] += 1
                        else:
                            fails.append(f'n={n} {code}: compat tree not isomorphic to T')
                    if labelled is not None:
                        Ap = [[sum(cof[cone_labels[S[i]]][j] * Binv[j][kk] for j in range(n))
                               for kk in range(n)] for i in range(n)]
                        # numeric sanity check of the conjugation at one point
                        x = [Q(3 * i + 1, i + 2) for i in range(n)]
                        Ax = [sum(A[i][j] * x[j] for j in range(n)) for i in range(n)]
                        z = [sum(B[i][j] * x[j] for j in range(n)) for i in range(n)]
                        lhs = [sum(B[i][j] * x[j] * Ax[j] for j in range(n)) for i in range(n)]
                        rhs = [z[i] * sum(Ap[i][j] * z[j] for j in range(n)) for i in range(n)]
                        if lhs != rhs:
                            fails.append(f'n={n} {code}: conjugated matrix wrong')
                        for Tp in labelled:
                            cnt['labelled_target_trees_tested'] += 1
                            member = in_family(Ap, n, Tp)
                            allcompat = all(tuple(sorted(e)) in cg for e in Tp)
                            if member != allcompat:
                                fails.append(f'n={n} {code}: membership != compatibility for {Tp}')
                            if member:
                                cnt['family_memberships'] += 1
                                if V.tree_code(n, Tp) != code:
                                    fails.append(f'n={n} {code}: member of L_T\' with T\' not iso T')
        per_n[str(n)] = cnt
    return result(EXACT, fails, statement=(
        'For every cone spanning tree S at generic exact points (n<=6): the n-subsets of rays that are '
        'linearly independent are exactly the cone spanning trees; the compatibility graph on S equals '
        'the set of size-2 fundamental paths (computed linearly and checked against graph fundamental '
        'paths); whenever it is a spanning tree it is isomorphic to T. For n<=5, the S-conjugated matrix '
        'lies in L_{T\'} (labelled T\' via Pruefer codes) iff every edge of T\' is compatible, and then T\' is '
        'isomorphic to T.'), per_n=per_n)


def _components(n, es):
    parent = list(range(n))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for u, v in es:
        parent[find(u)] = find(v)
    return {find(v) for v in range(n)}


# --------------------------------------------------------------------------
# (c) Sharpness and strata
# --------------------------------------------------------------------------

def check_c():
    fails = []
    stats = {'factor_points': 0, 'zero_d_points': 0, 'equal_points': 0, 'opposite_points': 0,
             'pencil_points': 0}
    ray_counts = {}
    for n in range(2, 7):
        for code, edges in trees_upto(n)[n].items():
            adj = adjacency(n, edges)
            base_factors = [f for f, _ in delta_factors(n, edges, {e: 1 for e in oriented(edges)})]
            for fac in base_factors:
                for k in range(2):
                    rng = seeded('c', n, code, fac, k)
                    a, d = generic_point(n, edges, rng)
                    kind = fac[0]
                    if kind == 'd':
                        _, v, j = fac
                        d[(v, j)] = 0
                    elif kind == 'minus':
                        _, j, u, v = fac
                        d[(v, j)] = d[(u, j)]
                    else:
                        _, j, u, v = fac
                        d[(v, j)] = -d[(u, j)]
                    vals = delta_factors(n, edges, d)
                    zeros = [f for f, x in vals if x == 0]
                    if zeros != [fac]:
                        fails.append(f'n={n} {code} {fac}: construction zeros {zeros}')
                        continue
                    stats['factor_points'] += 1
                    A = tree_matrix(n, edges, a, d)
                    cr = complete_rays(A)
                    got = ray_set(cr)
                    if pencils(cr):
                        fails.append(f'n={n} {code} {fac}: unexpected pencil')
                    rays = tree_rays(n, edges, d)
                    if kind == 'd':
                        stats['zero_d_points'] += 1
                        _, v, j = fac
                        e0 = tuple(sorted((v, j)))
                        want = {(tuple(i for i in range(n) if vec[i] != 0), normalise(vec))
                                for lab, vec in rays.items()
                                if not (lab[0] == 'r' and tuple(sorted(lab[1:])) == e0)}
                        if got != want or len(got) != 2 * n - 2:
                            fails.append(f'n={n} {code} {fac}: zero-d stratum')
                    else:
                        _, j, u, v = fac
                        # orient P_uj = d_{j->u} x_u - d_{u->j} x_j
                        Puj = [Q(0)] * n
                        Puj[u], Puj[j] = Q(d[(j, u)]), Q(-d[(u, j)])
                        Pvj = [Q(0)] * n
                        Pvj[v], Pvj[j] = Q(d[(j, v)]), Q(-d[(v, j)])
                        if kind == 'minus':
                            stats['equal_points'] += 1
                            extra = [x - y for x, y in zip(Puj, Pvj)]
                            want_sup = tuple(sorted((u, v)))
                        else:
                            stats['opposite_points'] += 1
                            extra = list(Puj)
                            extra[v] -= d[(j, v)]
                            want_sup = tuple(sorted((u, j, v)))
                        sup = tuple(i for i in range(n) if extra[i] != 0)
                        want = expected_tree_rayset(n, edges, d) | {(sup, normalise(extra))}
                        if sup != want_sup or got != want or len(got) != 2 * n:
                            fails.append(f'n={n} {code} {fac}: {kind} stratum')
                    ray_counts.setdefault(f'{kind} n={n}', set()).add(len(got))
            # pencil: d_{u->v} = d_{v->u} = 0
            for u, v in edges:
                a, d = generic_point(n, edges, seeded('c-pencil', n, code, u, v))
                d[(u, v)] = d[(v, u)] = 0
                cr = complete_rays(A := tree_matrix(n, edges, a, d))
                pen = [(r['support'], r['dim']) for r in pencils(cr)]
                stats['pencil_points'] += 1
                if (tuple(sorted((u, v))), 2) not in pen:
                    fails.append(f'n={n} {code} edge {u}{v}: no pencil on the edge')
    # the iff statement on small-integer points
    iff = {}
    for n in range(2, 6):
        for code, edges in trees_upto(n)[n].items():
            rng = seeded('c-iff', n, code)
            cnt = {'points': 2000, 'delta_zero': 0, 'delta_nonzero': 0, 'tree_form': 0,
                   'mismatches': 0}
            for _ in range(2000):
                a = [rng.randint(-2, 2) for _ in range(n)]
                b = {e: rng.randint(-2, 2) for e in oriented(edges)}
                d = {(v, j): b[(v, j)] - a[j] for v, j in oriented(edges)}
                delta = delta_value(n, edges, d)
                tf = is_tree_form(n, edges, complete_rays(tree_matrix(n, edges, a, d)))
                cnt['delta_zero' if delta == 0 else 'delta_nonzero'] += 1
                cnt['tree_form'] += tf
                if tf != (delta != 0):
                    cnt['mismatches'] += 1
                    if cnt['mismatches'] <= 3:
                        fails.append(f'iff n={n} {code}: a={a} b={sorted(b.items())}')
            iff[f'n={n} {code}'] = cnt
    return result(EXACT, fails, statement=(
        'For every tree n<=6 and every factor of Delta_T, at exact points where exactly that factor '
        'vanishes: d_{v->j}=0 gives exactly the coordinates and the other n-2 edge rays (2n-2 rays); '
        'd_{u->j}=d_{v->j} gives D_T plus P_uj-P_vj with support {u,v} (2n rays); '
        'd_{u->j}=-d_{v->j} gives D_T plus a ray proportional to P_uj - d_{j->v} x_v with support {u,j,v} '
        '(2n rays); no pencils. d_{u->v}=d_{v->u}=0 gives a pencil on {u,v}. On 2000 small-integer points '
        'per tree (n<=5, entries in -2..2) tree form holds iff Delta_T != 0.'),
        stats=stats, ray_counts={k: sorted(v) for k, v in sorted(ray_counts.items())},
        iff_points=iff)


# --------------------------------------------------------------------------
# (d) Family intersections
# --------------------------------------------------------------------------

def family_intersection(n, edges, j, u, v):
    Tp = tuple(sorted([tuple(sorted(e)) for e in edges if set(e) != {u, j}] + [tuple(sorted((u, v)))]))
    is_tree = len(Tp) == n - 1 and len(_components(n, Tp)) == 1
    eT, eTp = lt_equations(n, edges), lt_equations(n, Tp)
    extra = [0] * (n * n)
    extra[u * n + j] = 1
    extra[v * n + j] = -1
    r_int = rank(eT + eTp, n * n)
    r_ext = rank(eT + [extra], n * n)
    r_all = rank(eT + eTp + [extra], n * n)
    return Tp, is_tree, n * n - rank(eT, n * n), n * n - r_int, r_int == r_ext == r_all


def check_d():
    fails = []
    n = 4
    edges = ((0, 1), (1, 2), (2, 3))
    Tp, is_tree, dimT, dimI, same = family_intersection(n, edges, 1, 0, 2)
    if Tp != ((0, 2), (1, 2), (2, 3)) or not is_tree or dimI != 9 or dimT != 10 or not same:
        fails.append('P4 example')
    a = [1, 2, 3, 4]
    b = {(1, 0): 5, (0, 1): 6, (2, 1): 6, (1, 2): 7, (3, 2): 8, (2, 3): 9}
    d = {(x, y): b[(x, y)] - a[y] for x, y in b}
    M = tree_matrix(n, edges, a, d)
    if not (in_family(M, n, edges) and in_family(M, n, Tp)) or all(x == 0 for r in M for x in r):
        fails.append('P4 integer point')
    general = {}
    cases = 0
    for n in range(3, 8):
        cnt = 0
        for code, edges in trees_upto(n)[n].items():
            adj = adjacency(n, edges)
            for j in range(n):
                for u, v in permutations(adj[j], 2):
                    Tp2, is_tree2, dT, dI, same2 = family_intersection(n, edges, j, u, v)
                    cnt += 1
                    if not (is_tree2 and dT == 3 * n - 2 and dI == 3 * n - 3 and same2):
                        fails.append(f'n={n} {code} j={j} u={u} v={v}')
        general[str(n)] = cnt
        cases += cnt
    return result(EXACT, fails, statement=(
        'P4 = 0-1-2-3 and T\' = {02,12,23}: L_T cap L_T\' = {A in L_T : b_{0->1} = b_{2->1}}, '
        'dim L_T = 10, dim of intersection 9 = 3n-3. For all trees n<=7, every vertex j and ordered pair '
        'of distinct neighbours (u,v): T\' = T - uj + uv is a tree, dim(L_T cap L_T\') = 3n-3 and the '
        'intersection is cut out in L_T by b_{u->j} = b_{v->j}.'),
        p4_example={'T': [[0, 1], [1, 2], [2, 3]], 'T_prime': [list(e) for e in Tp],
                    'dim_L_T': dimT, 'dim_intersection': dimI,
                    'integer_point_params': {'a': a, 'b': {f'b_{x}->{y}': b[(x, y)] for x, y in sorted(b)}},
                    'integer_point_matrix': [vstr(r) for r in M]},
        general_cases_per_n=general, general_cases_total=cases)


# --------------------------------------------------------------------------
# (e) Normalised family = tangent space
# --------------------------------------------------------------------------

def check_e():
    fails = []
    per_tree = []
    for n in range(2, 7):
        monos = [(p, q) for p in range(n) for q in range(p, n)]
        midx = {m: k for k, m in enumerate(monos)}
        N = n * len(monos)

        def var(k, m):
            return k * len(monos) + midx[m]
        for code, edges in trees_upto(n)[n].items():
            rows = []
            for i in range(n):
                for m in monos:
                    if i not in m:
                        row = [0] * N
                        row[var(i, m)] = 1
                        rows.append(row)
            for u, v in edges:
                images = {}
                for m in monos:
                    img = tuple(sorted(v if t == u else t for t in m))
                    images.setdefault(img, []).append(m)
                for img in sorted(images):
                    row = [0] * N
                    for m in images[img]:
                        row[var(u, m)] += 1
                        row[var(v, m)] -= 1
                    rows.append(row)
            tangent = nullspace(rows, N)
            # LV normalised family {A in L_T : d_{u->v} = d_{v->u}}
            eqs = lt_equations(n, edges)
            for u, v in edges:
                row = [0] * (n * n)
                row[u * n + v] += 1
                row[v * n + v] -= 1
                row[v * n + u] -= 1
                row[u * n + u] += 1
                eqs.append(row)
            fam = nullspace(eqs, n * n)
            lv_vecs = []
            for vec in fam:
                X = [Q(0)] * N
                for k in range(n):
                    for j in range(n):
                        X[var(k, tuple(sorted((k, j))))] += vec[k * n + j]
                lv_vecs.append(X)
            rt, rl, ru = rank(tangent, N), rank(lv_vecs, N), rank(tangent + lv_vecs, N)
            if not (rt == rl == ru == 2 * n - 1):
                fails.append(f'n={n} {code}: dims {rt},{rl},{ru}')
            per_tree.append({'n': n, 'tree_code': code, 'unknowns': N,
                             'tangent_dim': rt, 'normalised_LV_dim': rl, 'sum_dim': ru})
    return result(EXACT, fails, statement=(
        'For every tree n<=6, the homogeneous quadratic vector fields on K^n tangent to every hyperplane '
        'x_i=0 and x_u-x_v=0 (uv in E) form a space of dimension 2n-1 equal to the LV fields with '
        'A in L_T and d_{u->v}=d_{v->u} on every edge.'), trees=per_tree)


# --------------------------------------------------------------------------
# (f) Caterpillar formula
# --------------------------------------------------------------------------

def R(ws):
    k = len(ws)
    if k == 0:
        return 1
    M = [[0] * k for _ in range(k)]
    for i in range(k):
        M[i][i] = ws[i] + (i > 0) + (i < k - 1)
        if i + 1 < k:
            M[i][i + 1] = M[i + 1][i] = -1
    return V.exact_det(M)


def caterpillar_q(lams):
    m = len(lams)
    w = [1 + (x > 0) for x in lams]
    if lams != lams[::-1]:
        return R(w)
    if m % 2 == 0:
        fix = R(w[:m // 2])
    else:
        c = (m - 1) // 2
        fix = w[c] * sum(R(w[:j]) for j in range((m - 1) // 2 + 1))
    total = R(w) + fix
    assert total % 2 == 0
    return total // 2


def spine(T):
    n = len(T)
    core = [v for v in range(n) if len(T[v]) >= 2]
    cs = set(core)
    if len(core) < 2:
        return None
    cdeg = {v: sum(w in cs for w in T[v]) for v in core}
    if max(cdeg.values()) > 2:
        return None
    start = min(v for v in core if cdeg[v] == 1)
    order, prev = [start], None
    while True:
        nxt = [w for w in T[order[-1]] if w in cs and w != prev]
        if not nxt:
            break
        prev = order[-1]
        order.append(nxt[0])
    return [sum(len(T[w]) == 1 for w in T[v]) for v in order]


def check_f():
    fails = []
    per_n = {}
    for n, Ts in C.trees_up_to(10):
        if n < 2:
            continue
        cnt = {'caterpillars': 0, 'core_matches': 0, 'orbit_checked': 0, 'orbit_matches': 0,
               'non_caterpillars': 0}
        for T in Ts:
            lams = spine(T)
            if lams is None:
                cnt['non_caterpillars'] += 1
                continue
            cnt['caterpillars'] += 1
            qf = caterpillar_q(lams)
            qc, _ = C.core_count(T)
            if qf == qc:
                cnt['core_matches'] += 1
            else:
                fails.append(f'n={n} lambda={lams}: formula {qf} core {qc}')
            if n <= 9:
                es = tuple(sorted(C.edges(T)))
                sizes, _ = V.orbit_counts(n, es, V.bases(n, es))
                cnt['orbit_checked'] += 1
                if len(sizes) == qf:
                    cnt['orbit_matches'] += 1
                else:
                    fails.append(f'n={n} lambda={lams}: formula {qf} orbits {len(sizes)}')
        per_n[str(n)] = cnt
    return result(EXACT, fails, statement=(
        'For every caterpillar (core a path v_1..v_m, m>=2) with n<=10, the closed formula q '
        '(R(w) = det(L_{P_k}+diag w), with the palindrome correction) equals the decorated-core '
        'count of lv_core_count.py, and for n<=9 equals the cone-basis orbit count of verify.py.'),
        per_n=per_n,
        caterpillars_total=sum(v['caterpillars'] for v in per_n.values()))


# --------------------------------------------------------------------------
# (g) Worked P3 example
# --------------------------------------------------------------------------

def check_g():
    import sympy as sp
    fails = []
    a1, a2, a3, b21, b12, b32, b23 = sp.symbols('a1 a2 a3 b_21 b_12 b_32 b_23')
    x = sp.symbols('x1:4')
    A = sp.Matrix([[a1, b12, b23], [b21, a2, b23], [b21, b32, a3]])
    d21, d12, d32, d23 = b21 - a1, b12 - a2, b32 - a2, b23 - a3
    P12 = d21 * x[0] - d12 * x[1]
    P23 = d32 * x[1] - d23 * x[2]
    delta = sp.factor(d21 * d12 * d32 * d23 * (d12 - d32) * (d12 + d32))
    Ax = A * sp.Matrix(x)
    C12 = a1 * x[0] + a2 * x[1] + b23 * x[2]
    C23 = b21 * x[0] + a2 * x[1] + a3 * x[2]
    for P, Cf, nm in ((P12, C12, 'P12'), (P23, C23, 'P23')):
        Pdot = sum(sp.diff(P, x[i]) * x[i] * Ax[i] for i in range(3))
        if sp.expand(Pdot - P * Cf) != 0:
            fails.append(f'symbolic Darboux identity {nm}')
    # concrete point
    an = [0, 0, 0]
    bn = {(1, 0): 1, (0, 1): 2, (2, 1): 3, (1, 2): 4}   # zero-based (v, j) for b_{v->j}
    dn = {e: bn[e] - an[e[1]] for e in bn}
    edges = ((0, 1), (1, 2))
    M = tree_matrix(3, edges, an, dn)
    dval = delta_value(3, edges, dn)
    cr = complete_rays(M)
    rays = sorted(ray_set(cr))
    tr = tree_rays(3, edges, dn)
    if [vstr(tr[('r', 0, 1)]), vstr(tr[('r', 1, 2)])] != [['1', '-2', '0'], ['0', '3', '-4']]:
        fails.append('integer edge rays')
    if len(rays) != 5 or pencils(cr) or dval == 0:
        fails.append('concrete point rays')
    names = {}
    for sup, vec in rays:
        if len(sup) == 1:
            names[(sup, vec)] = f'x{sup[0] + 1}'
        else:
            names[(sup, vec)] = 'P' + ''.join(str(s + 1) for s in sup)
    lab_paper = {('s', 0): 'x1', ('s', 1): 'x2', ('s', 2): 'x3', ('r', 0, 1): 'P12', ('r', 1, 2): 'P23'}
    trays = tree_rays(3, edges, dn)
    cone = V.cone_edges(3, edges)
    cone_labels = [('s', 0), ('s', 1), ('s', 2), ('r', 0, 1), ('r', 1, 2)]
    bases = V.bases(3, edges)
    tables = [p for p in permutations(range(3))]
    basis_records = []
    classes = {}
    for mask in sorted(bases):
        S = [j for j in range(5) if mask >> j & 1]
        hyper = fundamental_hypergraph(3, edges, mask)
        canon = min(tuple(sorted(sum(1 << p[s] for s in range(3) if h >> s & 1) for h in hyper))
                    for p in tables)
        classes.setdefault(canon, []).append(mask)
        Bm = [trays[cone_labels[j]] for j in S]
        Binv = inverse(Bm)
        z_supports = []
        for g in range(5):
            if mask >> g & 1:
                continue
            coeffs = [sum(trays[cone_labels[g]][t] * Binv[t][s] for t in range(3)) for s in range(3)]
            z_supports.append({'ray': lab_paper[cone_labels[g]],
                               'in_z_coordinates': vstr(coeffs)})
        basis_records.append({
            'basis_rays': [lab_paper[cone_labels[j]] for j in S],
            'cone_edges_paper_labels': [
                f'{{0,{cone[j][0] + 1}}}' if cone[j][1] == 3 else f'{{{cone[j][0] + 1},{cone[j][1] + 1}}}'
                for j in S],
            'fundamental_supports_z_indices': [[s + 1 for s in range(3) if h >> s & 1] for h in hyper],
            'nonbasis_rays_in_z_coordinates': z_supports})
    class_list = []
    for canon, masks in sorted(classes.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        class_list.append({'canonical_support_masks': list(canon),
                           'bases': [[lab_paper[cone_labels[j]] for j in range(5) if m >> j & 1]
                                     for m in masks]})
    if len(bases) != 8 or len(classes) != 2:
        fails.append('P3 bases/classes')
    rays_out = []
    for sup, vec in rays:
        rec = next(r for r in cr if r['support'] == sup)
        rays_out.append({'support_paper_labels': [s + 1 for s in sup],
                         'coefficients_x1_x2_x3': vstr(vec), 'cofactor_x1_x2_x3': vstr(rec['cofactor'])})
    return result(EXACT, fails, statement=(
        'Worked P3 example (paper labels 1-2-3, zero-based 0-1-2): symbolic matrix, edge Darboux '
        'polynomials and Delta; concrete integer point with 5 rays, 8 cone bases, 2 hypergraph classes.'),
        symbolic={'A': [[str(e) for e in A.row(i)] for i in range(3)],
                  'P12': str(sp.expand(P12)), 'P23': str(sp.expand(P23)),
                  'cofactor_P12': str(C12), 'cofactor_P23': str(C23),
                  'Delta': str(delta),
                  'Delta_in_d': 'd_{2->1} d_{1->2} d_{3->2} d_{2->3} (d_{1->2} - d_{3->2}) (d_{1->2} + d_{3->2})',
                  'd': {'d_{2->1}': 'b_21 - a1', 'd_{1->2}': 'b_12 - a2', 'd_{3->2}': 'b_32 - a2',
                        'd_{2->3}': 'b_23 - a3'}},
        concrete={'params': {'a1': 0, 'a2': 0, 'a3': 0, 'b_2->1': 1, 'b_1->2': 2, 'b_3->2': 3,
                             'b_2->3': 4},
                  'd': {'d_2->1': 1, 'd_1->2': 2, 'd_3->2': 3, 'd_2->3': 4},
                  'A': [vstr(r) for r in M], 'Delta': dval,
                  'edge_rays_integer_form': {'P12': 'x1 - 2*x2', 'P23': '3*x2 - 4*x3'},
                  'rays': rays_out, 'ray_count': len(rays),
                  'cone_bases': basis_records, 'cone_basis_count': len(bases),
                  'hypergraph_classes': class_list, 'hypergraph_class_count': len(classes)})


# --------------------------------------------------------------------------
# (h) van der Kamp's 13-parameter system (exploratory)
# --------------------------------------------------------------------------

def vdk_matrix(p):
    a = [None] + [Q(x) for x in p]
    return [[a[1], a[2], a[11], a[12], a[3]],
            [a[4], a[5], a[11], a[12], a[3]],
            [a[6], a[2], a[7], a[12], a[3]],
            [a[4], a[8], a[11], a[9], a[3]],
            [a[1] - (a[1] - a[4]) * (a[5] - a[10]) / (a[2] - a[5]), a[10], a[11], a[12], a[13]]]


def rank_function(vecs):
    m = len(vecs)
    return {mask: rank([vecs[i] for i in range(m) if mask >> i & 1]) if mask else 0
            for mask in range(1 << m)}


def matroid_isomorphic(r1, r2, m):
    """Backtracking matroid isomorphism using the full rank functions."""
    order = list(range(m))
    phi = [None] * m
    used = [False] * m

    def ok(k):
        # check all subsets of assigned elements {0..k} containing k
        for sub in range(1 << k):
            mask1 = sub | (1 << k)
            mask2 = 0
            for i in range(k + 1):
                if mask1 >> i & 1:
                    mask2 |= 1 << phi[i]
            if r1[mask1] != r2[mask2]:
                return False
        return True

    def rec(k):
        if k == m:
            return True
        for t in range(m):
            if not used[t]:
                phi[k] = t
                used[t] = True
                if ok(k) and rec(k + 1):
                    return True
                used[t] = False
        return False
    return rec(0)


def check_h():
    cone_ranks = {}
    for code, edges in trees_upto(5)[5].items():
        d = {e: 1 for e in oriented(edges)}
        # generic exact point for the cone matroid representation
        a, d = generic_point(5, edges, seeded('h-cone', code))
        vecs = [vec for _, vec in sorted(tree_rays(5, edges, d).items())]
        cone_ranks[code] = (edges, rank_function(vecs))
    # self-test of the isomorphism routine: each cone matroid at a second generic point,
    # with elements shuffled, is isomorphic to itself and to no other cone matroid
    selftest = {}
    for code, edges in sorted(trees_upto(5)[5].items()):
        a, d = generic_point(5, edges, seeded('h-cone-2', code))
        vecs = [vec for _, vec in sorted(tree_rays(5, edges, d).items())]
        seeded('h-shuffle', code).shuffle(vecs)
        r2 = rank_function(vecs)
        selftest[code] = {c2: matroid_isomorphic(r2, cone_ranks[c2][1], 9) for c2 in sorted(cone_ranks)}
    selftest_ok = all(v[c2] == (c2 == c) for c, v in selftest.items() for c2 in v)
    points = []
    rng = seeded('h-vdk')
    while len(points) < 5:
        p = [rng.randint(-9, 9) for _ in range(13)]
        if p[1] == p[4]:   # a2 == a5
            continue
        points.append(p)
    records = []
    for p in points:
        A = vdk_matrix(p)
        cr = complete_rays(A)
        rays = sorted(ray_set(cr))
        pens = [{'support': [s + 1 for s in r['support']], 'dim': r['dim']} for r in pencils(cr)]
        vecs = [list(v) for _, v in rays]
        rk = rank(vecs) if vecs else 0
        circ3 = [[[s + 1 for s in rays[i][0]] for i in t] for t in combinations(range(len(vecs)), 3)
                 if rank([vecs[i] for i in t]) < 3]
        rec = {'params_a1_to_a13': p, 'matrix': [vstr(r) for r in A],
               'ray_count': len(rays), 'pencils': pens,
               'rays': [{'support': [s + 1 for s in sup], 'coefficients': vstr(v)} for sup, v in rays],
               'rank': rk, 'three_circuits': circ3, 'three_circuit_count': len(circ3)}
        iso = {}
        if len(rays) == 9 and rk == 5:
            r1 = rank_function(vecs)
            rec['circuit_size_profile'] = circuit_profile(r1, 9)
            for code, (edges, r2) in sorted(cone_ranks.items()):
                iso[code] = matroid_isomorphic(r1, r2, 9)
            rec['isomorphic_to_cone_matroid_of_tree'] = iso
        else:
            rec['isomorphic_to_cone_matroid_of_tree'] = 'not tested (needs 9 rays of rank 5)'
        records.append(rec)
    return {'label': EXPLORATORY, 'status': 'PASS' if selftest_ok else 'FAIL',
            'note': ('Exploratory: records the complete linear Darboux structure of matrix (34) of '
                     'arXiv:2411.18264v4 at 5 deterministic integer points; status PASS means the '
                     'computation completed, not that any expectation holds.'),
            'cone_tree_codes': sorted(cone_ranks),
            'cone_circuit_size_profiles': {c: circuit_profile(r, 9) for c, (_, r) in sorted(cone_ranks.items())},
            'isomorphism_routine_selftest': selftest,
            'isomorphism_routine_selftest_status': 'PASS' if selftest_ok else 'FAIL',
            'points': records}


def circuit_profile(r, m):
    prof = {}
    for mask in range(1, 1 << m):
        k = bin(mask).count('1')
        if r[mask] == k - 1 and all(r[mask & ~(1 << i)] == k - 1 for i in range(m) if mask >> i & 1):
            prof[str(k)] = prof.get(str(k), 0) + 1
    return prof


# --------------------------------------------------------------------------
# (i) Kahan identity, (j) constant growth rate, (k) orthant
# --------------------------------------------------------------------------

def check_i():
    import sympy as sp
    fails = []
    x = sp.symbols('x1:4')
    y = sp.symbols('y1:4')
    h = sp.Symbol('h')

    def B(A, u, w):
        Au, Aw = A * sp.Matrix(u), A * sp.Matrix(w)
        return [(u[i] * Aw[i] + w[i] * Au[i]) / 2 for i in range(3)]

    def lin(c, u):
        return sum(ci * ui for ci, ui in zip(c, u))

    def checks(A, c, cof, name):
        Px, Py = lin(c, x), lin(c, y)
        Cx, Cy = lin(cof, x), lin(cof, y)
        Ax = A * sp.Matrix(x)
        Q_ = sp.expand(lin(c, [x[i] * Ax[i] for i in range(3)]) - Px * Cx)
        if Q_ != 0:
            fails.append(f'{name}: not a Darboux polynomial')
        pol = sp.expand(lin(c, B(A, x, y)) - (Px * Cy + Py * Cx) / 2)
        if pol != 0:
            fails.append(f'{name}: polarisation')
        kahan = sp.expand(Py * (1 - h * Cx / 2) - Px * (1 + h * Cy / 2)
                          - (lin(c, [y[i] - x[i] for i in range(3)]) - h * lin(c, B(A, x, y))))
        if kahan != 0:
            fails.append(f'{name}: Kahan rearrangement')

    a1, a2, a3, b21, b12, b32, b23 = sp.symbols('a1 a2 a3 b_21 b_12 b_32 b_23')
    A = sp.Matrix([[a1, b12, b23], [b21, a2, b23], [b21, b32, a3]])
    rays = {'x1': ([1, 0, 0], list(A.row(0))), 'x2': ([0, 1, 0], list(A.row(1))),
            'x3': ([0, 0, 1], list(A.row(2))),
            'P12': ([b21 - a1, -(b12 - a2), 0], [a1, a2, b23]),
            'P23': ([0, b32 - a2, -(b23 - a3)], [b21, a2, a3])}
    for name in sorted(rays):
        checks(A, rays[name][0], rays[name][1], f'P3 {name}')
    # general 3x3 matrix with its coordinate Darboux polynomials
    m = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'm{i + 1}{j + 1}'))
    for i in range(3):
        c = [int(k == i) for k in range(3)]
        checks(m, c, list(m.row(i)), f'general x{i + 1}')
    # general linear P, C: polarisation defect equals polarisation of the Darboux defect
    p = sp.symbols('p1:4')
    q = sp.symbols('q1:4')
    mx, my = m * sp.Matrix(x), m * sp.Matrix(y)

    def defect(u):
        mu = m * sp.Matrix(u)
        return lin(p, [u[i] * mu[i] for i in range(3)]) - lin(p, u) * lin(q, u)
    s = [x[i] + y[i] for i in range(3)]
    lhs = lin(p, B(m, x, y)) - (lin(p, x) * lin(q, y) + lin(p, y) * lin(q, x)) / 2
    if sp.expand(lhs - (defect(s) - defect(x) - defect(y)) / 2) != 0:
        fails.append('general polarisation of the Darboux defect')
    return result(EXACT, fails, statement=(
        'Symbolically (sympy): for every ray P of the generic P3 system and for the coordinate rays of a '
        'general 3x3 matrix, P(B(x,y)) = (P(x)C(y)+P(y)C(x))/2 and '
        'P(y)(1-hC(x)/2) - P(x)(1+hC(y)/2) = P(y-x) - h P(B(x,y)) identically, so the Kahan relation '
        'y - x = h B(x,y) gives P(x\')(1-hC(x)/2) = P(x)(1+hC(x\')/2). For general linear P, C and a general '
        '3x3 matrix, P(B(x,y)) - (P(x)C(y)+P(y)C(x))/2 is the polarisation of P(x*Ax) - P(x)C(x).'))


def check_j():
    fails = []
    edges = ((0, 1), (1, 2))
    a, d = generic_point(3, edges, seeded('j'))
    A = tree_matrix(3, edges, a, d)
    r0, r7 = complete_rays(A), complete_rays(A, r=[7, 7, 7])
    s0 = sorted((rr['support'], rr.get('vector'), rr['dim']) for rr in r0)
    s7 = sorted((rr['support'], rr.get('vector'), rr['dim']) for rr in r7)
    if s0 != s7:
        fails.append('ray sets differ')
    for rr in r7:
        if rr['dim'] == 1 and not darboux_identity_holds(A, rr['vector'], rr['cofactor'], r=[7, 7, 7], C0=7):
            fails.append(f'identity fails for {rr["support"]}')
    return result(EXACT, fails, statement=(
        'For a generic exact P3 point, the complete rays of x\' = x(7 + Ax) equal those of x\' = x(Ax) '
        '(cofactor shifted by the constant 7), each verified coefficientwise.'),
        ray_count=len(ray_set(r7)), params={'a': a, 'd': {f'd_{v}->{j}': d[(v, j)] for v, j in sorted(d)}})


def check_k():
    fails = []
    n, edges = 4, ((0, 1), (1, 2), (2, 3))
    cone_labels = [('s', v) for v in range(n)] + [('r', u, v) for u, v in edges]
    out = {}
    for signed in (False, True):
        a, d = generic_point(n, edges, seeded('k', signed), signed=signed)
        rays = tree_rays(n, edges, d)
        passing = []
        bases = V.bases(n, edges)
        for mask in sorted(bases):
            S = [j for j in range(2 * n - 1) if mask >> j & 1]
            B = [rays[cone_labels[j]] for j in S]
            ok_any = False
            for signs in product((1, -1), repeat=n):
                DB = [[s * x for x in row] for s, row in zip(signs, B)]
                if any(x < 0 for row in DB for x in row):
                    continue
                inv = inverse(DB)
                if all(x >= 0 for row in inv for x in row):
                    ok_any = True
                    break
            monomial = all(sum(x != 0 for x in row) == 1 for row in B)
            if ok_any != monomial:
                fails.append(f'signed={signed} basis {S}: orthant test != monomial')
            if ok_any:
                passing.append([label_str(cone_labels[j]) for j in S])
        if passing != [['s0', 's1', 's2', 's3']]:
            fails.append(f'signed={signed}: passing bases {passing}')
        out['signed_d' if signed else 'positive_d'] = {
            'bases_tested': len(bases), 'passing_bases': passing,
            'd': {f'd_{v}->{j}': d[(v, j)] for v, j in sorted(d)}}
    return result(EXACT, fails, statement=(
        'For P4 at a generic positive and a generic mixed-sign exact point, among the 21 cone spanning-tree '
        'bases only the all-spokes basis admits row signs D with DB >= 0 and (DB)^{-1} >= 0 '
        '(equivalently, only it is monomial).'), cases=out)


CHECKS = [('a_cofactor_lemma_and_triangles', check_a),
          ('b_one_sided_rigidity', check_b),
          ('c_sharpness_and_strata', check_c),
          ('d_family_intersection', check_d),
          ('e_normalised_family_tangent_space', check_e),
          ('f_caterpillar_formula', check_f),
          ('g_worked_P3_example', check_g),
          ('h_van_der_Kamp_13_parameter_system', check_h),
          ('i_kahan_identity', check_i),
          ('j_constant_growth_rate', check_j),
          ('k_orthant_bases', check_k)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=str(PACKAGE_ROOT / 'data' / 'extensions_verification.json'))
    ap.add_argument('--only', default=None, help='comma-separated check prefixes (debugging)')
    args = ap.parse_args()
    checks = {}
    for name, fn in CHECKS:
        if args.only and not any(name.startswith(p) for p in args.only.split(',')):
            continue
        checks[name] = fn()
        print(f'{name}: {checks[name]["status"]}', flush=True)
    exact_ok = all(c['status'] == 'PASS' for c in checks.values())
    output = {'status': 'PASS' if exact_ok else 'FAIL',
              'arithmetic': 'Python integers and fractions.Fraction; sympy for symbolic identities',
              'seeds': 'random.Random with fixed string seeds',
              'checks': checks}
    Path(args.out).write_text(json.dumps(output, indent=2, sort_keys=True) + '\n')
    print(f'overall: {output["status"]}; wrote {Path(args.out).name}', flush=True)
    return 0 if exact_ok else 1


if __name__ == '__main__':
    sys.exit(main())
