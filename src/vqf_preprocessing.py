import sympy as sp
from math import ceil
import itertools


# =========================== Section II: clauses, preprocessing ===========================
def build_clauses(m, n_p=None, n_q=None):
    """
    [Eq. 1-3] Builds the variables and one clause per bit position of m.
    Sizes (paper default): n_p = n_m, n_q = ceil(n_m/2); pass n_p, n_q to use known lengths.
    z[(i, t)] = carry from bit i into bit t (t > i); clause i is
    sum_j q_j p_{i-j} + incoming carries - m_i - sum_t 2^(t-i) z_{i,t}  (= 0).
    Returns p bits, q bits, carry variables, and the list of clause expressions.
    """
    nm = m.bit_length()
    n_p, n_q = n_p or nm, n_q or ceil(nm / 2)
    n_c = max(n_p + n_q - 1, nm)     # = n_p+n_q-1 (paper); covers all bits of m if n_p, n_q are given
    p = [sp.Symbol(f"p{i}") for i in range(n_p)]
    q = [sp.Symbol(f"q{i}") for i in range(n_q)]
    z = {(i, t): sp.Symbol(f"z_{i}_{t}") for i in range(n_c) for t in range(i + 1, n_c)}
    C = []
    for i in range(n_c):
        c = sp.Integer(0)
        for j in range(i + 1):
            if j < n_q and i - j < n_p:
                c += q[j] * p[i - j]                                      # sum_j q_j p_{i-j}
        c += sum((z[(j, i)] for j in range(i)), sp.Integer(0))            # incoming carries
        c -= (m >> i) & 1                                                 # m_i
        c -= sum((2 ** (t - i) * z[(i, t)] for t in range(i + 1, n_c)), sp.Integer(0))
        C.append(c)
    return p, q, z, C

def poly(expr):
    """
    [Helper] sympy expression -> {frozenset(variables): integer coefficient}.
    Variables are binary, so x^k = x and a monomial is just its set of variables.
    """
    d = {}
    for mono, c in sp.expand(expr).as_coefficients_dict().items():
        k = frozenset(mono.free_symbols)
        d[k] = d.get(k, 0) + int(c)
    return {k: c for k, c in d.items() if c}

def pmul(a, b):
    """
    [Section II helper] Multiply two polynomials in binary variables.

    A polynomial is a dict {frozenset(variables): integer coefficient}, e.g.
    p1*q1 - 1  ->  {frozenset({p1, q1}): 1, frozenset(): -1}.
    Because every variable is 0/1, x*x = x, so multiplying two monomials is just the
    union of their variable sets.
    Used in two places: energy() squares each clause with it (Eq. 6), and preprocess()
    uses it to expand substitutions such as y = 1 - x inside a clause.
    """
    out = {}
    for (k1, c1), (k2, c2) in itertools.product(a.items(), b.items()):
        k = k1 | k2
        out[k] = out.get(k, 0) + c1 * c2
    return {k: c for k, c in out.items() if c}

def preprocess(m, n_p=None, n_q=None, substitute=False, nontrivial=False):
    """
    [Section II.2] Simplifies the clauses: carry truncation + the Eq. 5 rules, repeated
    until nothing changes.
    substitute=True : also applies x - y = 0 => y = x and x + y - 1 = 0 => y = 1 - x
                      (exact reformulation; replaces Eq. 5 rule 2).
    nontrivial=True : also enforces p != 1 and q != 1 (removes the trivial factorisation).
    Returns p, q, z, known (var -> 0/1), sub (var -> (rep, flip): var = rep XOR flip),
    and the remaining clauses.
    """
    p, q, z, C = build_clauses(m, n_p, n_q)
    clauses = [poly(c) for c in C]
    known, zero, sub = {}, set(), {}         # sub[v] = (rep, flip)  <=>  v = rep XOR flip
    if n_p and n_q:                          # footnote 1: known lengths => leading bits are 1
        known[p[-1]] = known[q[-1]] = 1

    def lit(v):
        if v in known:
            return {frozenset(): 1} if known[v] else {}
        r, f = sub.get(v, (v, False))
        return {frozenset([r]): -1, frozenset(): 1} if f else {frozenset([r]): 1}

    def reduce(cl):
        out = {}
        for mono, c in cl.items():
            t = {frozenset(): c}
            for v in mono:
                t = pmul(t, lit(v))
            for k, cc in t.items():
                out[k] = out.get(k, 0) + cc
        return {k: c for k, c in out.items() if c and not any(zm <= k for zm in zero)}

    while True:
        before = (dict(known), set(zero), dict(sub))
        for v, (r, f) in list(sub.items()):                       # representative became known
            if r in known:
                known[v] = known[r] ^ int(f)
                del sub[v]
        for zm in list(zero):                                     # refresh xy=0 facts
            zero.discard(zm)
            if any(known.get(v) == 0 for v in zm):
                continue
            zm = frozenset(v for v in zm if known.get(v) != 1)
            if not zm:
                raise ValueError("inconsistent")
            if len(zm) == 1:
                known[next(iter(zm))] = 0
            else:
                zero.add(zm)
        if nontrivial:                                            # p != 1 and q != 1
            for bits in (p, q):
                hi = bits[1:]
                if any(known.get(v) == 1 for v in hi):
                    continue
                L = {sub.get(v, (v, False)) for v in hi if v not in known}
                if not L:
                    raise ValueError("only the trivial factorisation")
                if len(L) == 1:
                    (r, f), = L
                    known[r] = 1 - int(f)                         # literal r XOR f must be 1
        clauses = [r for r in (reduce(cl) for cl in clauses) if r]

        for cl in clauses:
            c0 = cl.get(frozenset(), 0)
            nz = {k: c for k, c in cl.items() if k}

            for k, c in nz.items():                               # carry truncation
                if len(k) == 1 and c < 0:
                    rest_max = c0 + sum(cc for kk, cc in nz.items() if kk != k and cc > 0)
                    if -c > rest_max:
                        known[next(iter(k))] = 0
            if not nz:
                raise ValueError("inconsistent")

            if substitute and len(nz) == 2 and all(len(k) == 1 for k in nz):
                (k1, w1), (k2, w2) = nz.items()                   # x - y = 0  or  x + y - 1 = 0
                if (w1 == -w2 and c0 == 0) or (w1 == w2 and c0 == -w1):
                    x, y = sorted((next(iter(k1)), next(iter(k2))), key=str)
                    flip = w1 == w2                               # y = x   or   y = 1 - x
                    sub[y] = (x, flip)
                    for v, (r, f) in list(sub.items()):
                        if r == y:
                            sub[v] = (x, f ^ flip)
                    break

            cs = list(nz.values())                                # Eq. (5) rules
            if all(c > 0 for c in cs) or all(c < 0 for c in cs):
                s = 1 if cs[0] > 0 else -1
                a, w = s * c0, [s * c for c in cs]
                if a > 0 or a < -sum(w):
                    raise ValueError("inconsistent")
                single = all(len(k) == 1 for k in nz)
                if a == 0 and single:                                         # rule 4
                    for k in nz:
                        known[next(iter(k))] = 0
                elif len(set(w)) == 1 and a == -w[0] * len(w):
                    if len(nz) == 1 and len(next(iter(nz))) == 2 or single:  # rules 1, 3, 5
                        for k in nz:
                            for v in k:
                                known[v] = 1
                elif single and len(w) == 2 and len(set(w)) == 1 and a == -w[0]:   # rule 2
                    zero.add(frozenset().union(*nz))
        if (known, zero, sub) == before:
            break

    clauses = [r for r in (reduce(cl) for cl in clauses) if r]
    return p, q, z, known, sub, clauses
