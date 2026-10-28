import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def load_results(filepath):
    with open(filepath) as f:
        return json.load(f)

# ── Plot 1: Runtime comparison (all three methods) ────────────────────────────
def plot_runtime(classical_results, shor_results, vqf_sim_results, vqf_hw_results=None):
    fig, ax = plt.subplots(figsize=(10, 5))

    for results, label, marker in [
        (classical_results, 'Classical', 'o'),
        (shor_results,      'Shor',      's'),
        (vqf_sim_results,   'VQF Sim',   '^'),
    ]:
        Ns = [r['N'] for r in results]
        ts = [r['runtime_total_sec'] for r in results]
        ax.plot(Ns, ts, marker=marker, label=label)

    if vqf_hw_results:
        Ns = [r['N'] for r in vqf_hw_results]
        ts = [r['runtime_total_sec'] for r in vqf_hw_results]
        ax.plot(Ns, ts, marker='D', label='VQF Hardware')

    ax.set_xlabel('N')
    ax.set_ylabel('Runtime (seconds)')
    ax.set_title('Runtime Comparison')
    ax.legend()
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig('figures/runtime_comparison.png', dpi=150)
    plt.show()

# ── Plot 2: Circuit depth and qubit count (Shor + VQF) ───────────────────────
def plot_circuit_metrics(shor_results, vqf_sim_results, vqf_hw_results=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for results, label, marker in [
        (shor_results,    'Shor',    's'),
        (vqf_sim_results, 'VQF Sim', '^'),
    ]:
        Ns     = [r['N'] for r in results if r['circuit_depth'] is not None]
        depths = [r['circuit_depth'] for r in results if r['circuit_depth'] is not None]
        qubits = [r['qubit_count'] for r in results if r['qubit_count'] is not None]
        ax1.plot(Ns, depths, marker=marker, label=label)
        ax2.plot(Ns, qubits, marker=marker, label=label)

    if vqf_hw_results:
        Ns     = [r['N'] for r in vqf_hw_results if r['circuit_depth'] is not None]
        depths = [r['circuit_depth'] for r in vqf_hw_results if r['circuit_depth'] is not None]
        qubits = [r['qubit_count'] for r in vqf_hw_results if r['qubit_count'] is not None]
        ax1.plot(Ns, depths, marker='D', label='VQF Hardware')
        ax2.plot(Ns, qubits, marker='D', label='VQF Hardware')

    ax1.set_xlabel('N'); ax1.set_ylabel('Circuit Depth'); ax1.set_title('Circuit Depth vs N')
    ax2.set_xlabel('N'); ax2.set_ylabel('Qubit Count');   ax2.set_title('Qubit Count vs N')
    ax1.legend(); ax2.legend()
    plt.tight_layout()
    plt.savefig('figures/circuit_metrics.png', dpi=150)
    plt.show()

# ── Plot 3: VQF success probability vs N ─────────────────────────────────────
def plot_success_prob(vqf_sim_results, vqf_hw_results=None):
    fig, ax = plt.subplots(figsize=(8, 5))

    Ns    = [r['N'] for r in vqf_sim_results]
    probs = [r['success_prob'] for r in vqf_sim_results]
    ax.plot(Ns, probs, marker='^', label='VQF Sim')

    if vqf_hw_results:
        Ns    = [r['N'] for r in vqf_hw_results]
        probs = [r['success_prob'] for r in vqf_hw_results]
        ax.plot(Ns, probs, marker='D', label='VQF Hardware')

    ax.set_xlabel('N')
    ax.set_ylabel('Success Probability')
    ax.set_title('VQF Success Probability vs N')
    ax.set_ylim(0, 1)
    ax.legend()
    plt.tight_layout()
    plt.savefig('figures/success_probability.png', dpi=150)
    plt.show()

# ── Plot 4: Success probability vs QAOA layers s (N=15) ──────────────────────
def plot_success_vs_layers(vqf_layer_sweep_results):
    """
    vqf_layer_sweep_results: list of Result dicts from running
    factor_vqf(15, cfg) with cfg.s = 1,2,3,4
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    layers = [r['layers_s'] for r in vqf_layer_sweep_results]
    probs  = [r['success_prob'] for r in vqf_layer_sweep_results]
    E_opts = [r['E_opt'] for r in vqf_layer_sweep_results]

    ax.plot(layers, probs, marker='o', label='Success Probability')
    ax2 = ax.twinx()
    ax2.plot(layers, E_opts, marker='s', linestyle='--', color='orange', label='E_opt')

    ax.set_xlabel('QAOA Layers (s)')
    ax.set_ylabel('Success Probability')
    ax2.set_ylabel('E_opt')
    ax.set_title('Success Probability vs QAOA Layers (N=15)')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('figures/success_vs_layers.png', dpi=150)
    plt.show()

# ── Plot 5: Gate count comparison (Shor + VQF hardware) ──────────────────────
def plot_gate_count(shor_results, vqf_hw_results):
    fig, ax = plt.subplots(figsize=(8, 5))

    for results, label, marker in [
        (shor_results,   'Shor',         's'),
        (vqf_hw_results, 'VQF Hardware', 'D'),
    ]:
        Ns     = [r['N'] for r in results if r['gate_count'] is not None]
        counts = [r['gate_count'] for r in results if r['gate_count'] is not None]
        ax.plot(Ns, counts, marker=marker, label=label)

    ax.set_xlabel('N')
    ax.set_ylabel('Gate Count')
    ax.set_title('Gate Count vs N')
    ax.legend()
    plt.tight_layout()
    plt.savefig('figures/gate_count.png', dpi=150)
    plt.show()
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def load_results(filepath):
    with open(filepath) as f:
        return json.load(f)

# ── Plot 1: Runtime comparison (all three methods) ────────────────────────────
def plot_runtime(classical_results, shor_results, vqf_sim_results, vqf_hw_results=None):
    fig, ax = plt.subplots(figsize=(10, 5))

    for results, label, marker in [
        (classical_results, 'Classical', 'o'),
        (shor_results,      'Shor',      's'),
        (vqf_sim_results,   'VQF Sim',   '^'),
    ]:
        Ns = [r['N'] for r in results]
        ts = [r['runtime_total_sec'] for r in results]
        ax.plot(Ns, ts, marker=marker, label=label)

    if vqf_hw_results:
        Ns = [r['N'] for r in vqf_hw_results]
        ts = [r['runtime_total_sec'] for r in vqf_hw_results]
        ax.plot(Ns, ts, marker='D', label='VQF Hardware')

    ax.set_xlabel('N')
    ax.set_ylabel('Runtime (seconds)')
    ax.set_title('Runtime Comparison')
    ax.legend()
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig('figures/runtime_comparison.png', dpi=150)
    plt.show()

# ── Plot 2: Circuit depth and qubit count (Shor + VQF) ───────────────────────
def plot_circuit_metrics(shor_results, vqf_sim_results, vqf_hw_results=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    for results, label, marker in [
        (shor_results,    'Shor',    's'),
        (vqf_sim_results, 'VQF Sim', '^'),
    ]:
        Ns     = [r['N'] for r in results if r['circuit_depth'] is not None]
        depths = [r['circuit_depth'] for r in results if r['circuit_depth'] is not None]
        qubits = [r['qubit_count'] for r in results if r['qubit_count'] is not None]
        ax1.plot(Ns, depths, marker=marker, label=label)
        ax2.plot(Ns, qubits, marker=marker, label=label)

    if vqf_hw_results:
        Ns     = [r['N'] for r in vqf_hw_results if r['circuit_depth'] is not None]
        depths = [r['circuit_depth'] for r in vqf_hw_results if r['circuit_depth'] is not None]
        qubits = [r['qubit_count'] for r in vqf_hw_results if r['qubit_count'] is not None]
        ax1.plot(Ns, depths, marker='D', label='VQF Hardware')
        ax2.plot(Ns, qubits, marker='D', label='VQF Hardware')

    ax1.set_xlabel('N'); ax1.set_ylabel('Circuit Depth'); ax1.set_title('Circuit Depth vs N')
    ax2.set_xlabel('N'); ax2.set_ylabel('Qubit Count');   ax2.set_title('Qubit Count vs N')
    ax1.legend(); ax2.legend()
    plt.tight_layout()
    plt.savefig('figures/circuit_metrics.png', dpi=150)
    plt.show()

# ── Plot 3: VQF success probability vs N ─────────────────────────────────────
def plot_success_prob(vqf_sim_results, vqf_hw_results=None):
    fig, ax = plt.subplots(figsize=(8, 5))

    Ns    = [r['N'] for r in vqf_sim_results]
    probs = [r['success_prob'] for r in vqf_sim_results]
    ax.plot(Ns, probs, marker='^', label='VQF Sim')

    if vqf_hw_results:
        Ns    = [r['N'] for r in vqf_hw_results]
        probs = [r['success_prob'] for r in vqf_hw_results]
        ax.plot(Ns, probs, marker='D', label='VQF Hardware')

    ax.set_xlabel('N')
    ax.set_ylabel('Success Probability')
    ax.set_title('VQF Success Probability vs N')
    ax.set_ylim(0, 1)
    ax.legend()
    plt.tight_layout()
    plt.savefig('figures/success_probability.png', dpi=150)
    plt.show()

# ── Plot 4: Success probability vs QAOA layers s (N=15) ──────────────────────
def plot_success_vs_layers(vqf_layer_sweep_results):
    """
    vqf_layer_sweep_results: list of Result dicts from running
    factor_vqf(15, cfg) with cfg.s = 1,2,3,4
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    layers = [r['layers_s'] for r in vqf_layer_sweep_results]
    probs  = [r['success_prob'] for r in vqf_layer_sweep_results]
    E_opts = [r['E_opt'] for r in vqf_layer_sweep_results]

    ax.plot(layers, probs, marker='o', label='Success Probability')
    ax2 = ax.twinx()
    ax2.plot(layers, E_opts, marker='s', linestyle='--', color='orange', label='E_opt')

    ax.set_xlabel('QAOA Layers (s)')
    ax.set_ylabel('Success Probability')
    ax2.set_ylabel('E_opt')
    ax.set_title('Success Probability vs QAOA Layers (N=15)')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig('figures/success_vs_layers.png', dpi=150)
    plt.show()

# ── Plot 5: Gate count comparison (Shor + VQF hardware) ──────────────────────
def plot_gate_count(shor_results, vqf_hw_results):
    fig, ax = plt.subplots(figsize=(8, 5))

    for results, label, marker in [
        (shor_results,   'Shor',         's'),
        (vqf_hw_results, 'VQF Hardware', 'D'),
    ]:
        Ns     = [r['N'] for r in results if r['gate_count'] is not None]
        counts = [r['gate_count'] for r in results if r['gate_count'] is not None]
        ax.plot(Ns, counts, marker=marker, label=label)

    ax.set_xlabel('N')
    ax.set_ylabel('Gate Count')
    ax.set_title('Gate Count vs N')
    ax.legend()
    plt.tight_layout()
    plt.savefig('figures/gate_count.png', dpi=150)
    plt.show()
