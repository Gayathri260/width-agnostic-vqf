from qiskit.quantum_info import SparsePauliOp
import itertools

def qubo_to_ising(E, qubits):
    """
    {frozenset(vars): coeff} -> SparsePauliOp via x = (1 - Z)/2, any degree.
    A monomial c * x_{v1}...x_{vd} expands to c/2^d * sum over subsets T of (-1)^|T| Z_T.
    Variable i is qubit i (rightmost-i character), matching bitstring[-(i+1)] in objective/decode.
    """
    n = len(qubits)
    idx = {v: i for i, v in enumerate(qubits)}
    terms = {}
    for mono, c in E.items():
        vs = sorted(idx[v] for v in mono)
        w = float(c) / 2 ** len(vs)
        for r in range(len(vs) + 1):
            for T in itertools.combinations(vs, r):
                terms[T] = terms.get(T, 0.0) + w * (-1) ** r
    out = []
    for T, c in terms.items():
        if abs(c) < 1e-12:
            continue
        lab = ['I'] * n
        for i in T:
            lab[n - 1 - i] = 'Z'
        out.append((''.join(lab), c))
    return SparsePauliOp.from_list(out) if out else SparsePauliOp('I' * n, coeffs=[0.0])
