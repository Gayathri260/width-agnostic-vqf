import json
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from collections import defaultdict

def load_results(filepath):
    with open(filepath) as f:
        return json.load(f)

# ── Plot 1: Runtime comparison by taking average of all run instances for an N (all three methods) ────────────────────────────
def plot_runtime(classical_results, shor_results, vqf_sim_results, vqf_hw_results=None):
    fig, ax = plt.subplots(figsize=(10, 5))

    def process_results(results):
        if not results:
            return [], []

        # Group runtimes by N
        grouped = defaultdict(list)
        for r in results:
            n = r.get('N')
            t = r.get('runtime_total_sec')
            if n is not None and t is not None:
                grouped[n].append(t)

        # Sort by N and calculate the average runtime for each N
        Ns = sorted(grouped.keys())
        avg_ts = [np.mean(grouped[n]) for n in Ns]
        return Ns, avg_ts

    datasets = [
        (classical_results, 'Classical', 'o'),
        (shor_results,      'Shor',      's'),
        (vqf_sim_results,   'VQF Sim',   '^'),
    ]

    for results, label, marker in datasets:
        Ns, ts = process_results(results)
        if Ns:
            ax.plot(Ns, ts, marker=marker, label=label)

    if vqf_hw_results:
        Ns, ts = process_results(vqf_hw_results)
        if Ns:
            ax.plot(Ns, ts, marker='D', label='VQF Hardware')

    ax.set_xlabel('N')
    ax.set_ylabel('Runtime (seconds)')
    ax.set_title('Runtime Comparison (Averaged by N)')
    ax.legend()
    ax.set_yscale('log')
    plt.tight_layout()
    plt.savefig('figures/runtime_comparison.png', dpi=150)
    plt.show()

# ── Plot 2: Circuit depth and qubit count (Shor + VQF) ───────────────────────
def plot_circuit_metrics(shor_results, vqf_sim_results, vqf_hw_results=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    width_unknown_vqf_sim = [
        r for r in vqf_sim_results
        if r.get('method') == "width-unknown-vqf-sim"
    ]
    for results, label, marker in [
        (shor_results,    'Shor',    's'),
        (width_unknown_vqf_sim, 'VQF Sim (Unknown Width)', '^'),
    ]:
        # Extract unique (N, depth, qubit) tuples to avoid duplicates/zigzags
        unique_runs = {}
        for r in results:
            n = r.get('N')
            depth = r.get('circuit_depth')
            qubit = r.get('qubit_count')
            if n is not None and depth is not None and qubit is not None:
                unique_runs[n] = (depth, qubit)

        if not unique_runs:
            continue
        sorted_ns = sorted(unique_runs.keys())
        depths = [unique_runs[n][0] for n in sorted_ns]
        qubits = [unique_runs[n][1] for n in sorted_ns]

        ax1.plot(sorted_ns, depths, marker=marker, label=label)
        ax2.plot(sorted_ns, qubits, marker=marker, label=label)

    if vqf_hw_results:
        Ns     = [r['N'] for r in vqf_hw_results if r['circuit_depth'] is not None]
        depths = [r['circuit_depth'] for r in vqf_hw_results if r['circuit_depth'] is not None]
        qubits = [r['qubit_count'] for r in vqf_hw_results if r['qubit_count'] is not None]
        ax1.plot(Ns, depths, marker='D', label='VQF Hardware')
        ax2.plot(Ns, qubits, marker='D', label='VQF Hardware')

    ax1.set_xlabel('N'); ax1.set_ylabel('Circuit Depth'); ax1.set_title('Circuit Depth vs N')
    ax1.set_yscale('log')
    ax2.set_xlabel('N'); ax2.set_ylabel('Qubit Count');   ax2.set_title('Qubit Count vs N')
    ax1.legend(); ax2.legend()
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax2.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig('figures/circuit_metrics.png', dpi=150)
    plt.show()

# ── Plot 3: VQF success probability vs N ─────────────────────────────────────
def plot_best_success_prob(best_results):
    """
    Plots the best success probabilities vs N for each category
    from the dictionary returned by get_best_success_probabilities().

    Args:
        best_results (dict): Nested dictionary structured as {category: {N: best_success_prob}}
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    # Define distinct visual styles for each category
    styles = {
        'width-unknown-vqf-sim': {
            'marker': '^',
            'linestyle': '-',
            'color': 'tab:blue',
            'label': 'Width-Unknown VQF Sim'
        },
        'width-known-vqf-sim': {
            'marker': 'o',
            'linestyle': '--',
            'color': 'tab:orange',
            'label': 'Width-Known VQF Sim'
        }
    }

    plotted_any = False

    for category, n_dict in best_results.items():
        if not n_dict:
            continue

        # Sort by N to ensure lines connect correctly across values
        Ns = sorted(n_dict.keys())
        probs = [n_dict[n] for n in Ns]

        # Fallback style if an unexpected category name appears
        style = styles.get(category, {
            'marker': 's',
            'linestyle': '-',
            'label': category
        })

        ax.plot(
            Ns, probs,
            marker=style['marker'],
            linestyle=style['linestyle'],
            color=style.get('color'),
            label=style['label']
        )
        plotted_any = True

    if not plotted_any:
        print("No valid data available to plot.")
        plt.close(fig)
        return

    ax.set_xlabel('N')
    ax.set_ylabel('Best Success Probability')
    ax.set_title('Best Success Probability vs N (Comparison)')
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig('figures/best_success_probability.png', dpi=150)
    plt.show()

# ── Plot 4: Best success probability vs QAOA layers s ──────────────────────
def plot_layer_sweep_best_success_prob(best_results_dict):
    """
    Plots the best success probability vs N for each layers_s value on a single graph.

    Args:
        best_results_dict (dict): Dictionary structured as {layers_s: {N: best_success_prob}}
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    # Styling options for different layer values (s = 1, 2, 3, 4)
    markers = ['o', 's', '^', 'D']
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

    plotted_any = False

    for idx, (layers, n_dict) in enumerate(sorted(best_results_dict.items())):
        if not n_dict:
            continue

        Ns = sorted(n_dict.keys())
        probs = [n_dict[n] for n in Ns]

        marker = markers[idx % len(markers)]
        color = colors[idx % len(colors)]

        ax.plot(
            Ns, probs,
            marker=marker,
            color=color,
            linestyle='-',
            label=f'{layers} layer(s)'
        )
        plotted_any = True

    if not plotted_any:
        print("No valid data available to plot.")
        plt.close(fig)
        return

    ax.set_xlabel('N')
    ax.set_ylabel('Best Success Probability')
    ax.set_title('Best Success Probability vs N across QAOA Layers')
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig('figures/layer_sweep_best_success_probability.png', dpi=150)
    plt.show()

# ── Plot 5: Runtime comparison per layers_s per N
def plot_layer_sweep_runtimes(results):
    """
    Calculates average runtimes per N for each layers_s value and plots them on a single graph.

    Args:
        results (list of dict): List of result dictionaries containing 'layers_s', 'N', and 'runtime_total_sec'.
    """
    fig, ax = plt.subplots(figsize=(9, 5))

    # Group runtimes by layers_s, then by N
    categorized_data = defaultdict(lambda: defaultdict(list))

    for r in results:
        n = r.get('N')
        layers = r.get('layers_s')
        runtime = r.get('runtime_total_sec')

        if n is None or layers is None or runtime is None:
            continue

        categorized_data[layers][n].append(runtime)

    markers = ['o', 's', '^', 'D']
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

    plotted_any = False

    for idx, (layers, n_dict) in enumerate(sorted(categorized_data.items())):
        if not n_dict:
            continue

        Ns = sorted(n_dict.keys())
        avg_runtimes = [float(np.mean(n_dict[n])) for n in Ns]

        marker = markers[idx % len(markers)]
        color = colors[idx % len(colors)]

        ax.plot(
            Ns, avg_runtimes,
            marker=marker,
            color=color,
            linestyle='-',
            label=f'Layers (s) = {layers}'
        )
        plotted_any = True

    if not plotted_any:
        print("No valid runtime data available to plot.")
        plt.close(fig)
        return

    ax.set_xlabel('N')
    ax.set_ylabel('Average Runtime (seconds)')
    ax.set_title('Runtime Comparison across QAOA Layers (Averaged by N)')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig('figures/layer_sweep_runtime_comparison.png', dpi=150)
    plt.show()

# ─ Plot 5: compare width-known VQF vs width-unknown VQF
def plot_comparison_best_prob(best_unknown_dict, best_known_dict):
    """
    Plots a comparison graph of the best success probabilities for
    width-unknown vs. width-known VQF simulation.

    Args:
        best_unknown_dict (dict): Dictionary of best probabilities for width-unknown {N: prob}
        best_known_dict (dict): Dictionary of best probabilities for width-known {N: prob}
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    plotted_any = False

    if best_unknown_dict:
        Ns_unknown = sorted(best_unknown_dict.keys())
        probs_unknown = [best_unknown_dict[n] for n in Ns_unknown]
        ax.plot(
            Ns_unknown, probs_unknown,
            marker='^', linestyle='-', color='tab:blue',
            label='Width-Unknown VQF Sim (Best)'
        )
        plotted_any = True

    if best_known_dict:
        Ns_known = sorted(best_known_dict.keys())
        probs_known = [best_known_dict[n] for n in Ns_known]
        ax.plot(
            Ns_known, probs_known,
            marker='o', linestyle='--', color='tab:orange',
            label='Width-Known VQF Sim (Best)'
        )
        plotted_any = True

    if not plotted_any:
        print("No data available to plot.")
        plt.close(fig)
        return

    ax.set_xlabel('N')
    ax.set_ylabel('Best Success Probability')
    ax.set_title('Comparison of Best Success Probabilities: Width-Unknown vs. Width-Known')
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig('figures/comparison_best_success_probability.png', dpi=150)
    plt.show()

# ── Plot 7: Gate count comparison (Shor + VQF hardware) ──────────────────────
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
