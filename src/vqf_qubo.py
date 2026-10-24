import sympy as sp

def energy(clauses):
    """[Eq. 6] E = sum of squared clauses, binary-reduced. Returns {frozenset(vars): coeff}."""
    E = {}
    for cl in clauses:
        for k, c in pmul(cl, cl).items():
            E[k] = E.get(k, 0) + c
    return {k: c for k, c in E.items() if c}

def reduce_to_quadratic(E):
    """
    For any monomial with degree > 2, introduce auxiliary variable
    w = x*y and add penalty 3w + xy - 2xw - 2yw
    Repeat until all terms are degree <= 2
    """
    changed = True
    while changed:
        changed = False
        new_E = {}
        aux_penalties = {}

        for mono, coeff in E.items():
            if len(mono) <= 2:
                new_E[mono] = new_E.get(mono, 0) + coeff
            else:
                # Pick first two vars to combine
                vars_list = sorted(mono, key=str)
                x, y = vars_list[0], vars_list[1]
                rest = frozenset(vars_list[2:])

                # Introduce w = x*y
                w = sp.Symbol(f'w_{x}_{y}')
                new_mono = frozenset([w]) | rest
                new_E[new_mono] = new_E.get(new_mono, 0) + coeff

                # Add penalty enforcing w = x*y (only once per pair)
                key = frozenset([x, y])
                if key not in aux_penalties:
                    P = 4  # penalty strength > max coefficient
                    aux_penalties[key] = w
                    # penalty: P*(3w + xy - 2xw - 2yw)
                    new_E[frozenset([w])] = new_E.get(frozenset([w]), 0) + 3*P
                    new_E[frozenset([x,y])] = new_E.get(frozenset([x,y]), 0) + P
                    new_E[frozenset([x,w])] = new_E.get(frozenset([x,w]), 0) - 2*P
                    new_E[frozenset([y,w])] = new_E.get(frozenset([y,w]), 0) - 2*P

                changed = True

        if changed:
            E = {k: v for k, v in new_E.items() if v != 0}

    return E

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
