import numpy as np
import scipy.linalg as la
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, Statevector, random_statevector, random_unitary
from verify_circuits import verify_unitary
import sys

def get_target_unitary(idx):
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    I = np.eye(2)
    
    if idx == 1:
        U = np.eye(4, dtype=complex)
        U[2:, 2:] = Y
        return U
    elif idx == 2:
        theta = np.pi/7
        Ry = la.expm(-1j * theta/2 * Y)
        U = np.eye(4, dtype=complex)
        U[2:, 2:] = Ry
        return U
    elif idx == 3:
        theta = np.pi/7
        return la.expm(1j * theta * np.kron(Z, Z))
    elif idx == 4:
        theta = np.pi/7
        H = np.kron(X, X) + np.kron(Y, Y)
        return la.expm(1j * theta * H)
    elif idx == 5:
        theta = np.pi/4
        H = np.kron(X, X) + np.kron(Y, Y) + np.kron(Z, Z)
        return la.expm(1j * theta * H)
    elif idx == 6:
        theta = np.pi/7
        H = np.kron(X, X) + np.kron(Z, I) + np.kron(I, Z)
        return la.expm(1j * theta * H)
    elif idx == 7:
        return None # Special case for state preparation
    elif idx == 8:
        return 0.5 * np.array([
            [1, 1, 1, 1],
            [1, 1j, -1, -1j],
            [1, -1, 1, -1],
            [1, -1j, -1, 1j]
        ])
    elif idx == 9:
        a = (1 - 1j) / 2
        b = (1 + 1j) / 2
        return np.array([
            [1, 0, 0, 0],
            [0, 0, -a, b],
            [0, 1j, 0, 0],
            [0, 0, -a, -b]
        ])
    elif idx == 10:
        return random_unitary(4, seed=42).data
    elif idx == 11:
        phases = [0, 4, 6, 5, 7, 7, 6, 5, 5, 7, 6, 6, 6, 6, 7, 5]
        diag = [np.exp(1j * p * np.pi/4) for p in phases]
        return np.diag(diag)
    return None

def verify_all():
    import os
    tasks = [
        (1, "task1_cy.qasm"),
        (2, "task2_cry.qasm"),
        (3, "task3_expzz.qasm"),
        (4, "task4_h1.qasm"),
        (5, "task5_h2.qasm"),
        (6, "task6_h3.qasm"),
        (7, "task7_state.qasm"),
        (8, "task8_u1.qasm"),
        (9, "task9_u2.qasm"),
        (10, "task10_random.qasm"),
        (11, "task11_diagonal.qasm")
    ]
    
    print(f"{'Task':<8} | {'T-count':<10} | {'Norm Distance':<15}")
    print("-" * 40)
    
    base_dir = os.path.join(os.path.dirname(__file__), "..", "circuits")
    for idx, filename in tasks:
        filepath = os.path.join(base_dir, filename)
        target = get_target_unitary(idx)
        if target is None:
            if idx == 7:
                psi_target = random_statevector(4, seed=42).data
                qc = QuantumCircuit.from_qasm_file(filepath)
                U = Operator(qc).data
                psi_actual = U[:, 0]
                fid = np.abs(np.vdot(psi_target, psi_actual))**2
                dist = np.sqrt(2 * (1 - np.sqrt(fid)))
                from verify_circuits import get_t_count
                t_count = get_t_count(qc)
                print(f"{idx:<8} | {t_count:<10} | {dist:.6e}")
            else:
                continue
        else:
            dist, t_count = verify_unitary(filepath, target, verbose=False)
            print(f"{idx:<8} | {t_count:<10} | {dist:.6e}")

if __name__ == "__main__":
    verify_all()
