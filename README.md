# Width-Agnostic VQF

Current implementations of VQF assume the bit-lengths of the factors the N. This is hardly applicable in real world, where RSA public keys are not distributed along with the factor bit-lengths. We have implemented a width-agnostic VQF algorithm, benchmarked on 6 odd, square-free bi-primes and compared with custom implementations of trial-division (classical) and Shor's simulation, all using the  Qiskit's quantum simulation framework.

## What is in here
- `src/classical.py`: trial division
- `src/shor.py`: Shor's algorithm
- `src/vqf_preprocessing.py`, `vqf_qubo.py`, `vqf_hamiltonian.py`, `vqf_sim.py`: VQF pipeline
- `src/config.py`: experiment configuration and result format; `src/utils.py`: utility functions; `src/plots.py`: plot functions
- `comparative_study.ipynb`: runs everything and produces the figures
- `results/`: JSON outputs; `figures/`: saved plots

## VQF pipeline
1. Binary long multiplication gives one equation per bit of N (Burges; Anschuetz et al.).
2. Classical preprocessing fixes or substitutes variables (Anschuetz Eq. 5; Karamlou et al.,
   App. C). Factor bit-lengths are not assumed known; only upper bounds derived from N are used.
3. Squared clauses give the cost function; Schaller-Schützhold penalty lowers the degree
   where possible.
4. Bits map to Pauli-Z operators, giving the Ising Hamiltonian.
5. QAOA with p layers; COBYLA minimizes ⟨H_C⟩, with Zhou et al. INTERP initialization.
6. Sampled bitstrings are decoded and checked; success is compared with random guessing.

## Setup
Note: This project was done solely with a research perspective and is only intended for academic research purposes.

```
git clone https://github.com/Gayathri260/quantum-factoring-comparative-study
cd quantum-factoring-comparative-study
pip install -r requirements.txt
jupyter notebook comparative_study.ipynb
```

Tested with Qiskit 2.4.2, Qiskit Aer 0.17.2, Python 3.14.7.

## References
- Anschuetz et al., Variational Quantum Factoring, QTOP 2019 (Springer); arXiv:1808.08927
- Karamlou et al., npj Quantum Inf. 7, 156 (2021); arXiv:2012.07825
- Hatwar, Aggarwal, Bhowmik, COMSNETS 2025; doi:10.1109/COMSNETS63942.2025.10885675
- Schaller & Schützhold, arXiv:0708.1882
- Xu et al., Phys. Rev. Lett. 108, 130501 (2012)
- Zhou et al., Phys. Rev. X 10, 021067 (2020)
- Farhi, Goldstone, Gutmann, arXiv:1411.4028
- Shor, Proc. FOCS 1994

## License
MIT
