import numpy as np
from collections import defaultdict

# prints summary of width-unknown runs
def get_best_success_prob(results):
    """
    Finds the best (maximum) success probability per N, separated by
    width-unknown-vqf-sim and width-known-vqf-sim based on bit_length.

    Args:
        results (list of dict): List of result dictionaries from the simulation.

    Returns:
        dict: Nested dictionary structured as {category: {N: best_success_prob}}
    """
    # Group success probabilities by category and N
    categorized_data = {
        'width-unknown-vqf-sim': defaultdict(list),
        'width-known-vqf-sim': defaultdict(list)
    }

    for r in results:
        prob = r.get('success_prob')
        n = r.get('N')

        # Skip entries without a success probability or N
        if prob is None or n is None:
            continue

        bit_len = r.get('bit_length')

        # Classify based on bit_length (width-known vs width-unknown)
        if bit_len == [None, None] or bit_len is None:
            category = 'width-unknown-vqf-sim'
        else:
            category = 'width-known-vqf-sim'

        categorized_data[category][n].append(prob)

    # Find the maximum (best) success probability per N for each category
    best_results = {}
    print("=== Best Success Probabilities per N ===")

    for category, n_dict in categorized_data.items():
        best_results[category] = {}
        print(f"\n[Category: {category}]")

        if not n_dict:
            print("  No entries found.")
            continue

        for n in sorted(n_dict.keys()):
            best_prob = max(n_dict[n])
            best_results[category][n] = best_prob
            print(f"  N = {n} -> {best_prob}")

    return best_results

def get_mean_success_prob(results):
    """
    Finds the mean (average) success probability per N, separated by
    width-unknown-vqf-sim and width-known-vqf-sim based on bit_length.

    Args:
        results (list of dict): List of result dictionaries from the simulation.

    Returns:
        dict: Nested dictionary structured as {category: {N: mean_success_prob}}
    """
    # Group success probabilities by category and N
    categorized_data = {
        'width-unknown-vqf-sim': defaultdict(list),
        'width-known-vqf-sim': defaultdict(list)
    }

    for r in results:
        prob = r.get('success_prob')
        n = r.get('N')

        # Skip entries without a success probability or N
        if prob is None or n is None:
            continue

        bit_len = r.get('bit_length')

        # Classify based on bit_length (width-known vs width-unknown)
        if bit_len == [None, None] or bit_len is None:
            category = 'width-unknown-vqf-sim'
        else:
            category = 'width-known-vqf-sim'

        categorized_data[category][n].append(prob)

    # Calculate the mean success probability per N for each category
    mean_results = {}
    print("=== Mean Success Probabilities per N ===")

    for category, n_dict in categorized_data.items():
        mean_results[category] = {}
        print(f"\n[Category: {category}]")

        if not n_dict:
            print("  No entries found.")
            continue

        for n in sorted(n_dict.keys()):
            probs = n_dict[n]
            mean_prob = float(np.mean(probs))
            mean_results[category][n] = mean_prob
            print(f"  N = {n} -> {mean_prob:.4f} (from {len(probs)} run(s))")

    return mean_results

def get_layer_sweep_mean_success_prob(results):
    """
    Calculates the mean success probability per N for each layers_s value.

    Args:
        results (list of dict): List of result dictionaries.

    Returns:
        dict: Nested dictionary structured as {layers_s: {N: mean_success_prob}}
    """
    categorized_data = defaultdict(lambda: defaultdict(list))

    for r in results:
        prob = r.get('success_prob')
        n = r.get('N')
        layers = r.get('layers_s')

        if prob is None or n is None or layers is None:
            continue

        categorized_data[layers][n].append(prob)

    mean_results = {}
    print("=== Layer Sweep: Mean Success Probabilities per N per QAOA Layer ===")

    for layers, n_dict in sorted(categorized_data.items()):
        mean_results[layers] = {}
        print(f"\n{layers} layer(s)")

        for n in sorted(n_dict.keys()):
            probs = n_dict[n]
            mean_prob = float(np.mean(probs))
            mean_results[layers][n] = mean_prob
            print(f"  N = {n} -> {mean_prob:.4f} (from {len(probs)} run(s))")

    return mean_results


def get_layer_sweep_best_success_prob(results):
    """
    Finds the best (maximum) success probability per N for each layers_s value.

    Args:
        results (list of dict): List of result dictionaries.

    Returns:
        dict: Nested dictionary structured as {layers_s: {N: best_success_prob}}
    """
    categorized_data = defaultdict(lambda: defaultdict(list))

    for r in results:
        prob = r.get('success_prob')
        n = r.get('N')
        layers = r.get('layers_s')

        if prob is None or n is None or layers is None:
            continue

        categorized_data[layers][n].append(prob)

    best_results = {}
    print("=== Layer Sweep: Best Success Probabilities per N per QAOA layer ===")

    for layers, n_dict in sorted(categorized_data.items()):
        best_results[layers] = {}
        print(f"\n{layers} layer(s)")

        for n in sorted(n_dict.keys()):
            best_prob = max(n_dict[n])
            best_results[layers][n] = best_prob
            print(f"  N = {n} -> {best_prob}")

    return best_results

from collections import defaultdict
import numpy as np

def get_width_known_mean_success_prob(results):
    """
    Calculates the mean success probability per N for width-known-vqf-sim.

    Args:
        results (list of dict): List of simulation result dictionaries.

    Returns:
        dict: Dictionary structured as {N: mean_success_prob}
    """
    n_dict = defaultdict(list)

    for r in results:
        method = r.get('method')
        bit_len = r.get('bit_length')
        is_width_known = (method == 'width-known-vqf-sim') or (bit_len is not None and bit_len != [None, None])

        if is_width_known:
            n = r.get('N')
            prob = r.get('success_prob')
            if n is not None and prob is not None:
                n_dict[n].append(prob)

    mean_results = {}
    print("=== Width-Known VQF Sim: Mean Success Probabilities per N ===")
    if not n_dict:
        print("No width-known-vqf-sim entries found.")

    for n in sorted(n_dict.keys()):
        probs = n_dict[n]
        mean_val = float(np.mean(probs))
        mean_results[n] = mean_val
        print(f"  N = {n} -> {mean_val:.4f} (from {len(probs)} run(s))")

    return mean_results

def get_width_known_best_success_prob(results):
    """
    Finds the best (maximum) success probability per N for width-known-vqf-sim.

    Args:
        results (list of dict): List of simulation result dictionaries.

    Returns:
        dict: Dictionary structured as {N: best_success_prob}
    """
    n_dict = defaultdict(list)

    for r in results:
        method = r.get('method')
        bit_len = r.get('bit_length')
        is_width_known = (method == 'width-known-vqf-sim') or (bit_len is not None and bit_len != [None, None])

        if is_width_known:
            n = r.get('N')
            prob = r.get('success_prob')
            if n is not None and prob is not None:
                n_dict[n].append(prob)

    best_results = {}
    print("=== Width-Known VQF Sim: Best Success Probabilities per N ===")
    if not n_dict:
        print("No width-known-vqf-sim entries found.")

    for n in sorted(n_dict.keys()):
        probs = n_dict[n]
        best_val = float(max(probs))
        best_results[n] = best_val
        print(f"  N = {n} -> {best_val:.4f}")

    return best_results
