import random
from sympy import isprime
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QFTGate, UnitaryGate

#====================== 1. Classical pre-check =========================
def is_prime_power(n:int) -> [bool, int|None, int|None] :
    for k in range(2, int(np.log2(n)+1)):
        p = round(n**(1/k))
        if round(p)**k == n and isprime(p):
            return (True, p, k)
    return (False, None, None)

def trivial_factor_check(n:int) -> [bool, [int, int]|None, str|None] :
    if n<2:
        raise ValueError("N must be >= 2")
    if n%2==0:
        return True, (2, n//2), "even"
    if isprime(n):
        return True, (1, n), "N is prime, nothing to factor"
    pp, p, k = is_prime_power(n)
    if pp:
        return True, (p, n//p), f"prime power: {p}^{k}"
    return False, None, None

# === 2. Shor - Classical part ===
def pick_random_a(n:int)->int :
    return random.randint(2, n-1)
def gcd(a:int, n:int)->int|None :
    d = math.gcd(a, n)
    if d!= 1:
        return d
    return None

# ============= 3. Shor- Quantum part (Finding order) ===============
def c_amodn(N:int, a:int, power:int, n:int) -> UnitaryGate :
    dim = 2 ** n
    a_pow = pow(a, power, N)
    perm = list(range(dim))
    for x in range(N):
        perm[x] = (a_pow * x) % N
    U = np.zeros((dim, dim))
    for x in range(dim):
        U[perm[x], x] = 1
    gate = UnitaryGate(U, label=f"{a}^{power} mod {N}")
    c_gate = gate.control()

    return c_gate

def build_circuit(N:int, a:int, m:int)->QuantumCircuit :
    n = math.ceil(np.log2(N)) #work qubits (ket psi)

    qc = QuantumCircuit(m+n, m)

    precision_qubits = list(range(m))
    work_qubits = list(range(m, m+n))

    #Initialize
    qc.h(precision_qubits)
    qc.x(work_qubits[0])

    #Controlled modular exponentiation
    for j in range(m):
        power = 2**j
        controlled_gate = c_amodn(N, a, power, n)
        qc.append(controlled_gate, [j]+work_qubits)

    #Apply Inverse QFT, then measure
    qc.append(QFTGate(m).inverse(), precision_qubits)
    qc.measure(precision_qubits, range(m))

    return qc

def run_circuit(qc:QuantumCircuit, shots:int = 1024)->dict[str,int] :
    backend = AerSimulator()
    transpiled = transpile(qc, backend)
    result = backend.run(transpiled, shots=shots).result()
    return result.get_counts()

# =================== 4. Shor - Classical post-processing ====================
def phase_to_order(measured_value:int, m:int, N:int) -> int|None :
    phase = measured_value/(2**m)
    frac = Fraction(phase).limit_denominator(N)
    r = frac.denominator
    return r
def order_to_factor(a:int, r:int, N:int) -> tuple[int]|None :
    '''
    r must be even
    a^(r/2) != -1 (mod N)  ==>  a^(r/2) != N-1
    '''
    if r%2 != 0:
        return None
    val = pow(a, r//2, N)

    if val==N-1:
        return None

    f1 = math.gcd(val-1, N)
    if f1>1 and f1<N:
        return f1, N//f1

    f2 = math.gcd(val+1, N)
    if f2>1 and f2<N:
        return f2, N//f2

    return None

# ================================ 5. Final Shor compiled in one function ==================================
def shor_factor(N:int, m:int=None, max_attempts:int=10, shots:int=1024) :
    '''
    N ->            integer to be factorised
    m ->            no. of QPE precision registers (default: no. of bits used to represent N in binary)
    max_attempts -> max number of attempts to try before quitting
    shots ->        no. of shots to take for measuring

    Returns (Result, QuantumCircuit)
    '''
    t_total = time.time()

    handled, factor, reason = trivial_factor_check(N)
    if handled:
        return Result(
                N, 'shor',None,
                True, factor, None,
                time.time()-t_total, None,
                None, None, None,
                None, None,
                None, None,
                None
            ), None #{"factors":factor, "method":reason, "attempt":0}
    if m is None:
        n = math.ceil(np.log2(N))
        m = 2*n + 1 #QPE precision register (at the top)

    t_exec = time.time()

    for attempt in range(1, max_attempts+1):
        a = pick_random_a(N)
        while gcd(a, N) != 1:
            a = pick_random_a(N)

        qc = build_circuit(N, a, m)
        counts = run_circuit(qc)
        result = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)

        for bitstring, freq in result:
            measured = int(bitstring, 2)
            r = phase_to_order(measured, m, N)
            if r is None:
                continue
            factors = order_to_factor(a, r, N)
            if factors is not None:
                transpiled_qc = transpile(qc, AerSimulator())
                return Result (
                        N, 'shor',None,
                        True, factors, 1/attempt,
                        time.time()-t_total, time.time()-t_exec,
                        transpiled_qc.depth(), qc.num_qubits, sum(transpiled_qc.count_ops().values()),
                        a, r,
                        None, None,
                        None
                    ), transpiled_qc, qc
            #{"factors":factors,"method":"shor","attempt":attempt,"a":a,"order":r,'circuit':qc,"circuit_depth": qc.depth(),"qubit_count": qc.num_qubits,"circuit_depth_transpiled": transpile(qc, AerSimulator()).depth()}
    return Result(
            N, 'shor',None,
            False, None, None,
            time.time()-t_total, time.time()-t_exec,
            None, None, None,
            None, None,
            None, None,
            None
        ), None, None #{"factors":None, "method":"failed", "attempt":attempt}
