import sympy as sp
from vqf_preprocessing import pmul
from fractions import Fraction

def clause_penalty(cl):
    """
    Penalty for one clause: 0 iff the clause equals 0, >= 1 otherwise.
    Plain square (Eq. 6), or, for a product AB with coefficient +-1, the Schaller-Schutzhold
    form: with clause = AB + S, penalty = ((2A + 2B - 1 + 4S)^2 - 1) / 8  [Sec. VI, Eq. 43].
    The lowest (max degree, term count) form is kept.
    """
    best = pmul(cl, cl)
    best_score = (max((len(k) for k in best), default=0), len(best))
    for ab, c in cl.items():
        if len(ab) != 2 or abs(c) != 1:
            continue
        A, B = sorted(ab, key=str)
        T = {frozenset([A]): 2, frozenset([B]): 2, frozenset(): -1}       # 2A + 2B - 1 + 4S
        for k, cc in cl.items():
            if k != ab:
                T[k] = T.get(k, 0) + 4 * c * cc
        T = {k: v for k, v in T.items() if v}
        sq = pmul(T, T)
        sq[frozenset()] = sq.get(frozenset(), 0) - 1
        cand = {k: Fraction(v, 8) for k, v in sq.items() if v}
        score = (max((len(k) for k in cand), default=0), len(cand))
        if score < best_score:
            best, best_score = cand, score
    return best

def energy(clauses):
    """[Eq. 6] Sum of clause penalties; identical (or sign-flipped) clauses are counted once."""
    E, seen = {}, set()
    for cl in clauses:
        key = frozenset(cl.items())
        neg = frozenset((k, -c) for k, c in cl.items())
        if key in seen or neg in seen:
            continue
        seen.add(key)
        for k, c in clause_penalty(cl).items():
            E[k] = E.get(k, 0) + c
    return {k: c for k, c in E.items() if c}


def qubit_list(E, p, q, known, sub):
    """
    [End of Section II: decides which variables become qubits]

    Returns the sorted list of variables that are still unknown after preprocessing:
    every variable that appears in the energy E, plus any p/q bit that is unknown but
    appears in no clause (for example p1 for N=15 after substitution, where both
    5x3 and 15x1 are valid). Qubit i corresponds to bit i of a basis-state index.
    Its length n is the qubit count of the algorithm.
    """
    free_pq = {sub.get(v, (v, False))[0] for v in p + q if v not in known}
    return sorted(set().union(*E.keys()) | free_pq, key=str)

def decode(assign, p, q, known, sub):
    """
    [Post-processing, after Section III: from a measured bit string back to factors]

    assign: {qubit variable: 0/1}, one measured bit string.
    Returns (P, Q), the two integers built from the bits of p and q.
    Each bit is resolved in this order: a known value (fixed by preprocessing), else
    its substitution (bit = representative XOR flip), else the measured value
    (default 0 if the variable is not in assign).
    Carries are not needed here, only p and q.
    """
    def val(v):
        if v in known:
            return known[v]
        r, f = sub.get(v, (v, False))
        return known.get(r, assign.get(r, 0)) ^ int(f)
    return (sum(val(v) << i for i, v in enumerate(p)),
            sum(val(v) << i for i, v in enumerate(q)))
