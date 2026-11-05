from dataclasses import dataclass, field
from typing import Optional, Tuple
import json, os
from datetime import datetime

@dataclass
class ExperimentConfig:
    shots: int = 4096
    n_restarts: int = 3
    maxiter_sim: int = 500       # COBYLA iterations for simulation
    maxiter_hw: int = 150        # QNSPSA iterations for hardware
    s: int = 2                   # QAOA layers
    mode: str = 'sim'            # 'sim' or 'hardware'
    backend_name: Optional[str] = None   # e.g. 'ibm_brisbane'; None = Aer
    optimization_level: int = 1
    p:int=None #factor bit_length
    q:int=None #factor bit_length

@dataclass
class Result:
    # Identity
    N: int
    method: str                          # 'classical', 'shor', 'vqf-sim', 'vqf-hardware'

    # Outcome
    success: bool
    factors: Optional[Tuple[int, int]]   # (p, q) or None
    success_prob: Optional[float]        # None for classical (deterministic)

    # Timing
    runtime_total_sec: float             # full pipeline including preprocessing
    runtime_exec_sec: Optional[float]    # circuit/QAOA execution only; None for classical

    # Circuit metrics
    circuit_depth: Optional[int]
    qubit_count: Optional[int]
    gate_count: Optional[int]

    # Shor-specific
    a: Optional[int]
    order_r: Optional[int]

    # VQF-specific
    E_opt: Optional[float]
    layers_s: Optional[int]
    bit_length:Optional[Tuple(int,int)] #factor bit-length p,q if inputted, else None

    # Run context
    backend: Optional[str]
    job_id:Optional[str] = None          # for hardware submissions
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

def save_result(result: Result, filepath: str):
    import dataclasses
    entry = dataclasses.asdict(result)
    if os.path.exists(filepath):
        with open(filepath) as f:
            all_results = json.load(f)
    else:
        all_results = []
    all_results.append(entry)
    with open(filepath, 'w') as f:
        json.dump(all_results, f, indent=2)

def load_results(filepath: str):
    with open(filepath) as f:
        return json.load(f)
