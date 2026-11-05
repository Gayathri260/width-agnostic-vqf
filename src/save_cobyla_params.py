import json
import os
import numpy as np

def save_cobyla_params(N, best_params, best_energy, layers_s, filepath='results/cobyla_params.json'):
    """
    Save best COBYLA parameters found during simulation.
    These will be used as initial point for hardware runs.
    """
    if best_params is None: return
    entry = {
        'N': N,
        'best_params': best_params.tolist(),
        'layers_s': layers_s,
        'order': 'beta[0..s-1] then gamma[0..s-1]',
        'best_energy': float(best_energy)
    }
    if os.path.exists(filepath):
        with open(filepath) as f:
            all_params = json.load(f)
    else:
        all_params = []

    # overwrite if N already exists
    all_params = [e for e in all_params if e['N'] != N]
    all_params.append(entry)

    with open(filepath, 'w') as f:
        json.dump(all_params, f, indent=2)
    print(f"  Params saved for N={N}")


def load_cobyla_params(filepath='results/cobyla_params.json'):
    with open(filepath) as f:
        entries = json.load(f)
    return {e['N']: np.array(e['best_params']) for e in entries}
