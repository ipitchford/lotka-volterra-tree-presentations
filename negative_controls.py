#!/usr/bin/env python3
"""Exact semantic negative controls (standard library only).

Each check verifies a statement of the paper on a small exact instance and then
confirms that a plausible corruption of the same statement is rejected.
Exit status is nonzero if any true statement fails or any corruption passes.
"""
import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def lv_rhs(A, n):
    # dx_i/dt = x_i (A x)_i as a dict {(i, j): coeff} of quadratic monomials x_i x_j (i <= j)
    out = []
    for i in range(n):
        poly = {}
        for j in range(n):
            key = tuple(sorted((i, j)))
            poly[key] = poly.get(key, 0) + Fraction(A[i][j])
        out.append(poly)
    return out


def is_darboux(A, c):
    """True if P = sum c_i x_i satisfies dP/dt = P * C for some linear C (exact linear algebra)."""
    n = len(c)
    f = lv_rhs(A, n)
    dP = {}
    for i in range(n):
        for k, v in f[i].items():
            dP[k] = dP.get(k, 0) + Fraction(c[i]) * v
    # unknown cofactor coefficients C_0..C_{n-1}; P*C has coefficient at (i,j): c_i C_j + c_j C_i (i<j), c_i C_i (i=j)
    rows, rhs = [], []
    for i in range(n):
        for j in range(i, n):
            row = [Fraction(0)] * n
            if i == j:
                row[i] += Fraction(c[i])
            else:
                row[j] += Fraction(c[i]); row[i] += Fraction(c[j])
            rows.append(row); rhs.append(dP.get((i, j), Fraction(0)))
    # solve least-squares-free: Gaussian elimination for consistency
    M = [r + [b] for r, b in zip(rows, rhs)]
    piv_row = 0
    for col in range(n):
        p = next((r for r in range(piv_row, len(M)) if M[r][col] != 0), None)
        if p is None:
            continue
        M[piv_row], M[p] = M[p], M[piv_row]
        for r in range(len(M)):
            if r != piv_row and M[r][col] != 0:
                fac = M[r][col] / M[piv_row][col]
                M[r] = [a - fac * b for a, b in zip(M[r], M[piv_row])]
        piv_row += 1
    return all(any(x != 0 for x in r[:n]) or r[n] == 0 for r in M)


failures = []


def expect(label, value, wanted):
    ok = value == wanted
    print(f"{'ok  ' if ok else 'FAIL'} {label}: {value} (expected {wanted})")
    if not ok:
        failures.append(label)


# Example 1.1: path 1-2-3 with a_j = 0, b(2->1)=1, b(1->2)=2, b(3->2)=3, b(2->3)=4
A = [[0, 2, 4], [1, 0, 4], [1, 3, 0]]
expect('Example 1.1: x1 - 2x2 is Darboux', is_darboux(A, [1, -2, 0]), True)
expect('Example 1.1: 3x2 - 4x3 is Darboux', is_darboux(A, [0, 3, -4]), True)
expect('corruption x1 - 3x2 rejected', is_darboux(A, [1, -3, 0]), False)
expect('corruption x1 + x3 (non-edge support) rejected', is_darboux(A, [1, 0, 1]), False)

# Theorem 2.3: exactly 2n-1 rays at a generic point; search all small integer forms up to scaling
def rays(A, bound=4):
    n = len(A); found = set()
    rng = range(-bound, bound + 1)
    def rec(prefix):
        if len(prefix) == n:
            c = list(prefix)
            if any(c):
                g = next(x for x in c if x != 0)
                key = tuple(Fraction(x, g) for x in c)
                if key not in found and is_darboux(A, c):
                    found.add(key)
            return
        for v in rng:
            rec(prefix + [v])
    rec([])
    return found
expect('Theorem 2.3: P3 example has 2n-1 = 5 rays (small search)', len(rays(A)), 5)

# Remark 2.4 (sharpness): P2 with d(1->2)=0, d(2->1)=1 has no edge ray
A2 = [[0, 0], [1, 0]]
expect('Remark 2.4: degenerate P2 has only the 2 coordinate rays', len(rays(A2)), 2)

# Proposition 3.8: matrix in L_P4 and in L_star(centre 3)
M = [[1, 6, 7, 9], [5, 2, 7, 9], [5, 6, 3, 9], [5, 6, 8, 4]]
def in_family(M, edges):
    n = len(M)
    return all(M[u][j] == M[v][j] for (u, v) in edges for j in range(n) if j not in (u, v))
expect('Proposition 3.8: matrix in L_P4', in_family(M, [(0, 1), (1, 2), (2, 3)]), True)
expect('Proposition 3.8: matrix in L_K13 (centre 3)', in_family(M, [(0, 2), (1, 2), (3, 2)]), True)
M_bad = [row[:] for row in M]; M_bad[0][1] = 7
expect('corruption b(1->2) != b(3->2) leaves L_K13', in_family(M_bad, [(0, 2), (1, 2), (3, 2)]), False)

# Theorem 7.1: q(P_n) = (F_{2n-2} + F_{n-1})/2 against the committed core-forest counts
F = [0, 1]
while len(F) < 40:
    F.append(F[-1] + F[-2])
data = json.loads((ROOT / 'data' / 'lv_core_counts.json').read_text())
paths = {}
for block in data['by_n']:
    n = block['n']
    for t in block['trees']:
        if sorted(t['edges']) == [[i, i + 1] for i in range(n - 1)]:
            paths[n] = t['q']
formula = {n: (F[2 * n - 2] + F[n - 1]) // 2 for n in paths}
expect('Theorem 7.1: path counts n<=10 match data', paths == formula and len(paths) == 9, True)
wrong = {n: (F[2 * n - 2] + F[n]) // 2 for n in paths}
expect('corruption (F_{2n-2}+F_n)/2 rejected', paths == wrong, False)

print('negative_controls: ' + ('PASS' if not failures else 'FAIL ' + ', '.join(failures)))
sys.exit(1 if failures else 0)
