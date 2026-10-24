from qiskit.quantum_info import SparsePauliOp

def qubo_to_ising(E, qubits):
    """
    Convert QUBO energy dict {frozenset(vars): coeff} to SparsePauliOp.

    Binary variable substitution: x = (1 - Z) / 2
    - Constant term c        → c * I
    - Linear term c*x_i      → c/2 * I - c/2 * Z_i
    - Quadratic term c*x_i*x_j → c/4 * I - c/4 * Z_i - c/4 * Z_j + c/4 * Z_iZ_j
    """
    n = len(qubits)
    var_index = {v: i for i, v in enumerate(qubits)}
    pauli_dict = {}

    def add(pauli_str, coeff):
        pauli_dict[pauli_str] = pauli_dict.get(pauli_str, 0) + coeff

    I = 'I' * n

    for mono, coeff in E.items():
        coeff = float(coeff)
        degree = len(mono)

        if degree == 0:
            add(I, coeff)

        elif degree == 1:
            v = list(mono)[0]
            i = var_index[v]
            add(I, coeff / 2)
            z = list(I); z[i] = 'Z'
            add(''.join(z), -coeff / 2)

        elif degree == 2:
            v1, v2 = list(mono)
            i, j = var_index[v1], var_index[v2]
            add(I, coeff / 4)
            zi = list(I); zi[i] = 'Z'; add(''.join(zi), -coeff / 4)
            zj = list(I); zj[j] = 'Z'; add(''.join(zj), -coeff / 4)
            zij = list(I); zij[i] = 'Z'; zij[j] = 'Z'; add(''.join(zij), coeff / 4)

        else:
            raise ValueError(f"Degree {degree} term found — reduce_to_quadratic not applied")

    pauli_list = [(p, c) for p, c in pauli_dict.items() if abs(c) > 1e-10]
    return SparsePauliOp.from_list(pauli_list)
