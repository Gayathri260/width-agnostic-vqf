import time
from config import ExperimentConfig, Result, save_result
from vqf_preprocessing import preprocess
from vqf_qubo import energy, reduce_to_quadratic, qubit_list, decode
from vqf_hamiltonian import qubo_to_ising

def factor_vqf(N: int, cfg: ExperimentConfig) -> tuple:
    t_total = time.time()

    # --- Preprocessing ---
    p, q, z, known, sub, clauses = preprocess(N)
    E = energy(clauses)
    E2 = reduce_to_quadratic(E)
    qubits = qubit_list(E2, p, q, known, sub)
    n = len(qubits)

    print(f"N={N}: {n} qubits after preprocessing")

    if n == 0:
        P, Q = decode({}, p, q, known, sub)
        return Result(
            N, 'vqf-preprocessing',
            P*Q == N and P > 1 and Q > 1, (P, Q), None,
            time.time()-t_total, None,
            None, None, None,
            None, None,
            cfg.backend_name, None,
            None
        ), None

    # --- Hamiltonian + Circuit ---
    hamiltonian = qubo_to_ising(E2, qubits)

    from qiskit.circuit.library import QAOAAnsatz
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    qaoa_circ = QAOAAnsatz(hamiltonian, reps=cfg.s, flatten=True)
    qaoa_circ.measure_all()

    if cfg.mode == 'sim':
        from qiskit_aer import AerSimulator
        backend = AerSimulator(method='matrix_product_state')
        pm = generate_preset_pass_manager(
            optimization_level=cfg.optimization_level,
            basis_gates=['cx', 'u', 'id', 'reset'],
            coupling_map=None
        )
    else:
        from qiskit_ibm_runtime import QiskitRuntimeService
        service = QiskitRuntimeService()
        backend = service.backend(cfg.backend_name)
        pm = generate_preset_pass_manager(
            optimization_level=cfg.optimization_level,
            backend=backend
        )

    transpiled_qc = pm.run(qaoa_circ)
    print(f"  Circuit depth: {transpiled_qc.depth()}, params: {transpiled_qc.num_parameters}")

    # --- Optimization ---
    import numpy as np
    from scipy.optimize import minimize

    t_exec = time.time()

    best_params, best_energy = None, float('inf')

    if cfg.mode == 'sim':
        def objective(params):
            bound = transpiled_qc.assign_parameters(params)
            job = backend.run(bound, shots=cfg.shots)
            counts = job.result().get_counts()
            total = sum(counts.values())
            exp_val = 0.0
            for bitstring, count in counts.items():
                bits = {v: int(bitstring[-(i+1)]) for i, v in enumerate(qubits)}
                e = 0.0
                for mono, coeff in E2.items():
                    term = 1
                    for var in mono:
                        term *= bits.get(var, 0)
                    e += coeff * term
                exp_val += (count / total) * e
            return exp_val
        for i in range(cfg.n_restarts):
            rng = np.random.default_rng(i)
            theta0 = rng.uniform(0, 2*np.pi, transpiled_qc.num_parameters)
            opt = minimize(objective, theta0, method='COBYLA',
                          options={'maxiter': cfg.maxiter_sim, 'rhobeg': 0.5})
            print(f"  restart {i+1}/{cfg.n_restarts}: E={opt.fun:.4f}")
            if opt.fun < best_energy:
                best_energy = opt.fun
                best_params = opt.x
    else:
        from qiskit_ibm_runtime import Session
        from qiskit_ibm_runtime.executor_sampler import Sampler
        from qiskit_algorithms.optimizers import QNSPSA

        with Session(backend=backend) as session:
            sampler = Sampler(mode=session)
            sampler.options.default_shots = cfg.shots

            def hw_objective(params):
                bound = transpiled_qc.assign_parameters(params)
                job = sampler.run([bound])
                counts = job.result()[0].data.meas.get_counts()
                total = sum(counts.values())
                exp_val = 0.0
                for bitstring, count in counts.items():
                    bits = {v: int(bitstring[-(i+1)]) for i, v in enumerate(qubits)}
                    e = 0.0
                    for mono, coeff in E2.items():
                        term = 1
                        for var in mono:
                            term *= bits.get(var, 0)
                        e += coeff * term
                    exp_val += (count / total) * e
                return exp_val

            for i in range(cfg.n_restarts):
                rng = np.random.default_rng(i)
                theta0 = rng.uniform(0, 2*np.pi, transpiled_qc.num_parameters)
                opt = QNSPSA(maxiter=cfg.maxiter_hw)
                result = opt.minimize(fun=hw_objective, x0=theta0)
                print(f"  restart {i+1}/{cfg.n_restarts}: E={result.fun:.4f}")
                if result.fun < best_energy:
                    best_energy = result.fun
                    best_params = result.x

    t_exec_total = time.time() - t_exec

    # --- Final sample with best params ---
    bound_final = transpiled_qc.assign_parameters(best_params)

    if cfg.mode == 'sim':
        final_counts = backend.run(bound_final, shots=cfg.shots).result().get_counts()
    else:
        with Session(backend=backend) as session:
            sampler = Sampler(mode=session)
            job = sampler.run([bound_final])
            final_counts = job.result()[0].data.meas.get_counts()

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
        print(f"  {bs}: {cnt}/{total} → {P}×{Q} valid={P*Q==N and P>1 and Q>1}")

    # --- Result ---
    return Result(
        N, f'vqf-{cfg.mode}',
        best_factors is not None, best_factors, success_prob,
        time.time() - t_total, t_exec_total,
        transpiled_qc.depth(), n, sum(transpiled_qc.count_ops().values()),
        best_energy, cfg.s,
        cfg.backend_name, None,
        None
    ), transpiled_qc
