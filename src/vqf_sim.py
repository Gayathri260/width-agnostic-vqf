import time
import json
import os
import numpy as np
from scipy.optimize import minimize

from config import ExperimentConfig, Result, save_result
from vqf_preprocessing import preprocess
from vqf_qubo import energy, qubit_list, decode
from vqf_hamiltonian import qubo_to_ising

def interp_init(prev, p):
    """Zhou et al., PRX 10, 021067, Eq. (B1): linearly interpolate the optimum at depth p
    to a starting point at depth p+1. prev = [beta_0..beta_{p-1}, gamma_0..gamma_{p-1}]."""
    def one(v):
        v = np.concatenate(([0.0], v, [0.0]))
        return [(i - 1) / p * v[i - 1] + (p - i + 1) / p * v[i] for i in range(1, p + 2)]
    return np.array(one(prev[:p]) + one(prev[p:]))

def factor_vqf_sim(N: int, cfg: ExperimentConfig):
    """
    Full VQF simulation pipeline.
    Returns (Result, transpiled_circuit).
    """
    from qiskit.circuit.library import QAOAAnsatz
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    from qiskit_aer import AerSimulator

    t_total = time.time()

    # --- Preprocessing ---
    p, q, z, known, sub, clauses = preprocess(N, cfg.p, cfg.q)
    E = energy(clauses)
    qubits = qubit_list(E, p, q, known, sub)
    n = len(qubits)

    print(f"N={N}: {n} qubits after preprocessing")

    if not E:
        P, Q = decode({}, p, q, known, sub)
        return Result(
            N, 'vqf-preprocessing',
            P*Q == N and P > 1 and Q > 1, (P, Q), None,
            time.time()-t_total, None,
            None, None, None,
            None, None,
            None, None,(cfg.p, cfg.q),
            None, None
        ), None, None, None #Result(), transpiled circuit, circuit, best_params

    # --- Hamiltonian + Circuit ---
    hamiltonian = qubo_to_ising(E, qubits)
    qaoa_circ = QAOAAnsatz(hamiltonian, reps=cfg.s, flatten=True)
    qaoa_circ.measure_all()

    backend = AerSimulator(method='statevector')
    pm = generate_preset_pass_manager(
        optimization_level=cfg.optimization_level,
        basis_gates=['cx', 'u', 'id', 'reset'],
        coupling_map=None
    )
    transpiled_qc = pm.run(qaoa_circ)
    print(f"  depth={transpiled_qc.depth()}, "
          f"params={transpiled_qc.num_parameters}, "
          f"gates={sum(transpiled_qc.count_ops().values())}")

    # --- COBYLA Optimization ---

    def objective(params, qc):
        bound = qc.assign_parameters(params)
        counts = backend.run(bound, shots=cfg.shots).result().get_counts()
        total = sum(counts.values())
        exp_val = 0.0
        for bitstring, count in counts.items():
            bits = {v: int(bitstring[-(i+1)]) for i, v in enumerate(qubits)}
            e = 0.0
            for mono, coeff in E.items():
                term = 1
                for var in mono:
                    term *= bits.get(var, 0)
                e += coeff * term
            exp_val += (count / total) * e
        return exp_val

    def build(reps):
        c = QAOAAnsatz(hamiltonian, reps=reps, flatten=True)
        c.measure_all()
        return pm.run(c)

    t_exec = time.time()
    best_params, best_energy = None, float('inf')

    prev = None
    for p_ in range(1, cfg.s + 1):
        qc_p = build(p_)
        if p_ == 1:                          # two parameters only: random multi-start
            best = None
            for i in range(cfg.n_restarts):
                rng = np.random.default_rng(i)
                r = minimize(objective, rng.uniform(0, 2 * np.pi, 2), args=(qc_p,), method='COBYLA',
                             options={'maxiter': cfg.maxiter_sim, 'rhobeg': 0.5})
                if best is None or r.fun < best.fun:
                    best = r
        else:                                # deeper layers start from the interpolated optimum
            best = minimize(objective, interp_init(prev, p_ - 1), args=(qc_p,), method='COBYLA',
                            options={'maxiter': cfg.maxiter_sim, 'rhobeg': 0.2})
        prev = best.x
        print(f"  p={p_}: E={best.fun:.4f}")
    best_params, best_energy = prev, best.fun
    transpiled_qc = qc_p

    t_exec_total = time.time() - t_exec

    # --- Final Sample ---
    bound_final = transpiled_qc.assign_parameters(best_params)
    final_counts = backend.run(
        bound_final, shots=cfg.shots
    ).result().get_counts()

    # --- Decode ---
    total = sum(final_counts.values())
    success_count = 0
    best_factors = None

    for bitstring, count in final_counts.items():
        bits = {v: int(bitstring[-(i+1)]) for i, v in enumerate(qubits)}
        P, Q = decode(bits, p, q, known, sub)
        if P * Q == N and P > 1 and Q > 1:
            success_count += count
            if best_factors is None:
                best_factors = (P, Q)

    success_prob = success_count / total

    # --- Diagnostic ---
    print(f"\n── Top 5 bitstrings (N={N}) ──")
    for bs, cnt in sorted(final_counts.items(), key=lambda x: -x[1])[:5]:
        bits = {v: int(bs[-(i+1)]) for i, v in enumerate(qubits)}
        P, Q = decode(bits, p, q, known, sub)
        print(f"  {bs}: {cnt}/{total} → {P}×{Q} "
              f"valid={P*Q==N and P>1 and Q>1}")

    return Result(
        N, 'width-known-vqf-sim' if cfg.p else 'width-unknown-vqf-sim',
        best_factors is not None, best_factors, success_prob,
        time.time()-t_total, t_exec_total,
        transpiled_qc.depth(), n, sum(transpiled_qc.count_ops().values()),
        None, None,
        best_energy, cfg.s, (cfg.p,cfg.q),
        None, None
    ), transpiled_qc, qaoa_circ, best_params
