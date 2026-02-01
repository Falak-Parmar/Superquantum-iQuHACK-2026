import json
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Pauli, PauliList, Clifford
import qiskit.qasm2
import os

def solve_bonus():
    # Adjusted paths for reorganization
    base_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(base_dir, '../docs/Challenge 12.json')
    output_path = os.path.join(base_dir, '../circuits/bonus_solution.qasm')
    
    if not os.path.exists(json_path):
        print(f"Error: Could not find {json_path}")
        return

    with open(json_path, 'r') as f:
        data = json.load(f)

    n = data['n']
    terms = data['terms']

    pauli_strings = [t['pauli'] for t in terms]
    ks = [t['k'] for t in terms]

    # Basis finding
    all_sym = []
    for s in pauli_strings:
        p = Pauli(s)
        all_sym.append(np.concatenate([p.x, p.z]).astype(int))
    all_sym = np.array(all_sym)
    
    rows, cols = all_sym.shape
    A = all_sym.copy().astype(int)
    r = 0
    mapping = np.arange(rows)
    for c in range(cols):
        if r >= rows: break
        pivot = r + np.argmax(A[r:, c])
        if A[pivot, c] == 0: continue
        A[[r, pivot]] = A[[pivot, r]]
        mapping[[r, pivot]] = mapping[[pivot, r]]
        for i in range(rows):
            if i != r and A[i, c] == 1:
                A[i] = (A[i] + A[r]) % 2
        r += 1
    
    basis_indices = mapping[:r]
    print(f"Basis rank: {len(basis_indices)}")
    
    diag_qc = QuantumCircuit(n)
    cur_basis = PauliList([pauli_strings[idx] for idx in basis_indices])
    
    def apply_gate(gate_name, qubits):
        g_qc = QuantumCircuit(n)
        if gate_name == 'h':
            diag_qc.h(qubits[0])
            g_qc.h(qubits[0])
        elif gate_name == 'sdg':
            diag_qc.sdg(qubits[0])
            g_qc.sdg(qubits[0])
        elif gate_name == 'cx':
            diag_qc.cx(qubits[0], qubits[1])
            g_qc.cx(qubits[0], qubits[1])
        elif gate_name == 'swap':
            diag_qc.swap(qubits[0], qubits[1])
            g_qc.swap(qubits[0], qubits[1])
        
        nonlocal cur_basis
        cur_basis = cur_basis.evolve(Clifford(g_qc))

    for i in range(len(basis_indices)):
        # 1. Pivot
        target_p = -1
        target_q = -1
        for q_idx in range(i, n):
            for p_idx in range(i, len(basis_indices)):
                if cur_basis[p_idx].x[q_idx] or cur_basis[p_idx].z[q_idx]:
                    target_p = p_idx
                    target_q = q_idx
                    break
            if target_p != -1: break
        
        if target_p == -1: continue
        if target_p != i:
            tmp = cur_basis[i]
            cur_basis[i] = cur_basis[target_p]
            cur_basis[target_p] = tmp
        if target_q != i:
            apply_gate('swap', [i, target_q])

        # 2. Local transformation to Z at qubit i
        if cur_basis[i].x[i]:
            if cur_basis[i].z[i]: # Y
                apply_gate('sdg', [i])
            apply_gate('h', [i])
        
        # 3. Clean Row i at qubits k != i
        # CX(k, i) maps Z_i -> Z_k Z_i.
        # If cur_basis[i] has Z_k and Z_i, apply CX(k, i) results in (Z_k Z_i) Z_k = Z_i.
        for k in range(n):
            if k == i: continue
            if cur_basis[i].x[k]:
                apply_gate('h', [k])
                apply_gate('cx', [k, i]) # control=k, target=i
                apply_gate('h', [k])
            if cur_basis[i].z[k]:
                apply_gate('cx', [k, i]) # control=k, target=i

        # 4. Clean column i for other rows j != i
        # CX(i, j) maps Z_j -> Z_i Z_j.
        # If cur_basis[j] has Z_i and Z_j, apply CX(i, j) results in (Z_i Z_j) Z_i = Z_j.
        for j in range(len(basis_indices)):
            if j == i: continue
            if cur_basis[j].z[i]:
                apply_gate('cx', [i, j]) # control=i, target=j

    print("Diagonalization completed.")
    for i in range(len(basis_indices)):
        print(f"Basis {i} -> {cur_basis[i].to_label()}")

    C = Clifford(diag_qc)
    transformed_all = PauliList(pauli_strings).evolve(C)
    
    final_qc = QuantumCircuit(n)
    final_qc.append(diag_qc, range(n))

    for i, p_trans in enumerate(transformed_all):
        if p_trans.to_label().endswith('I'*n): continue
        sign = 1
        if p_trans.phase == 2: sign = -1
        
        label_clean = p_trans.to_label()
        for p in ['-', 'i', '+']:
            if label_clean.startswith(p):
                label_clean = label_clean[1:]
        
        z_indices = [idx for idx, char in enumerate(reversed(label_clean)) if char == 'Z']
        if not z_indices: continue
        
        k_eff = ks[i] * sign
        for j in range(len(z_indices) - 1):
            final_qc.cx(z_indices[j], z_indices[-1])
        final_qc.rz(k_eff * np.pi / 4, z_indices[-1])
        for j in range(len(z_indices) - 2, -1, -1):
            final_qc.cx(z_indices[j], z_indices[-1])

    final_qc.append(diag_qc.inverse(), range(n))
    final_qc = transpile(final_qc, basis_gates=['h', 't', 'tdg', 'cx'], optimization_level=3)
    
    with open(output_path, 'w') as f:
        f.write(qiskit.qasm2.dumps(final_qc))
    print(f"Solution saved to {output_path}. T-count: {final_qc.count_ops().get('t', 0) + final_qc.count_ops().get('tdg', 0)}")

if __name__ == "__main__":
    solve_bonus()
