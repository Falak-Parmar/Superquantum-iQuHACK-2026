import json
import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, Pauli
import qiskit.qasm2
import scipy.linalg as la
import os

def verify_bonus():
    print("Verifying Bonus Challenge...")
    # Adjusted paths for reorganization
    json_path = os.path.join(os.path.dirname(__file__), '../docs/Challenge 12.json')
    qasm_path = os.path.join(os.path.dirname(__file__), '../circuits/bonus_solution.qasm')
    
    if not os.path.exists(json_path):
        print(f"Error: Could not find {json_path}")
        return
    if not os.path.exists(qasm_path):
        print(f"Error: Could not find {qasm_path}")
        return

    with open(json_path, 'r') as f:
        data = json.load(f)
    
    n = data['n']
    terms = data['terms']
    
    # Reconstruct target unitary
    U_target = np.eye(2**n, dtype=complex)
    for term in terms:
        p_str = term['pauli']
        k = term['k']
        P_mat = Pauli(p_str).to_matrix()
        angle = k * np.pi / 8
        Rot = la.expm(-1j * angle * P_mat)
        U_target = Rot @ U_target

    # Load solution
    try:
        qc = qiskit.qasm2.load(qasm_path)
        U_sim = Operator(qc).data
    except Exception as e:
        print(f"Failed to load/simulate solution: {e}")
        return

    # Compare
    fid = np.abs(np.trace(U_sim.conj().T @ U_target)) / (2**n)
    print(f"Fidelity (up to global phase): {fid}")
    
    # Normalize global phase
    phase = np.angle(np.trace(U_sim.conj().T @ U_target))
    U_sim_adj = U_sim * np.exp(1j * phase)
    dist_adj = np.linalg.norm(U_target - U_sim_adj)
    print(f"Adjusted operator distance: {dist_adj}")

    if fid > 0.9999 or dist_adj < 1e-8:
        print("VERIFICATION SUCCESSFUL")
    else:
        print("VERIFICATION FAILED")

    # Check gate set
    allowed = {'h', 't', 'tdg', 'cx'}
    ops = qc.count_ops()
    print(f"Ops: {ops}")
    for op in ops:
        if op not in allowed:
            print(f"WARNING: Disallowed gate found: {op}")

if __name__ == "__main__":
    verify_bonus()
