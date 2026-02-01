# Superquantum iQuHACK 2026 Challenge

This repository contains the solution for the Superquantum iQuHACK 2026 Challenge. The goal of the challenge was to synthesize various quantum unitaries into Clifford+T sequences with high precision.

## Repository Structure

- `solver.py`: The main script used to synthesize the target unitaries into Clifford+T QASM files.
- `verify_all.py`: A wrapper script to verify all 11 tasks and report their performance metrics.
- `verify_circuits.py`: Core logic for calculating operator norm distance and T-counts.
- `task*.qasm`: Generated Clifford+T circuits for each of the 11 challenge tasks.
- `.gitignore`: Standard git ignore file to keep the repository clean.

## Getting Started

### Prerequisites

- Python 3.8+
- Qiskit
- Scipy
- Numpy
- mpmath
- pygridsynth (for 1-qubit rotation synthesis)

### Running the Solver

To regenerate the Clifford+T circuits, run:

```bash
python3 solver.py
```

### Verification

To verify the generated circuits and see the performance statistics, run:

```bash
python3 verify_all.py
```

## Performance Results

All tasks were verified using the operator norm distance metric and T-count calculation.

| Task | Description | Norm Distance | T-count |
| :--- | :--- | :--- | :--- |
| 1 | Controlled-Y | 0.0000e+00 | 4 |
| 2 | Controlled-Ry(π/7) | 6.8145e-13 | 504 |
| 3 | exp(i π/7 ZZ) | 1.1120e-12 | 2956 |
| 4 | exp(i π/7 (XX+YY)) | 8.1579e-13 | 2758 |
| 5 | exp(i π/4 H2) | 2.2204e-16 | 0 |
| 6 | exp(i π/7 H3) | 1.4597e-12 | 2540 |
| 7 | State Preparation | 6.7878e-07 | 2522 |
| 8 | Structured 1 | 1.9131e-12 | 2674 |
| 9 | Structured 2 | 1.6163e-12 | 3332 |
| 10 | Random Unitary | 2.4678e-12 | 5126 |
| 11 | 4-qubit Diagonal | 2.6346e-12 | 3660 |

## Methodology

1. **Decomposition**: 2-qubit unitaries are decomposed into CX and 1-qubit U3 gates using Qiskit's transpiler with level 3 optimization.
2. **Synthesis**: 1-qubit U3 gates are further decomposed into Rz-H-Rz-H-Rz sequences. The Rz rotations are synthesized into Clifford+T sequences using `pygridsynth` for near-optimal T-count given a target precision.
3. **Verification**: The resulting QASM files are loaded back into Qiskit operators and compared against the target matrices using the operator norm distance, accounting for global phase.
