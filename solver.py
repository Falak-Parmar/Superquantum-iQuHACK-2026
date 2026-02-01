import numpy as np
import scipy.linalg as la
from mpmath import mp
import pygridsynth
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator, Statevector
import os

mp.dps = 50

def get_clifford_t_qasm(theta_phi_lambda, qubit):
    """
    Synthesizes a 1-qubit unitary U(theta, phi, lambda) using pygridsynth.
    Decomposes U into Rz(phi + pi/2) H Rz(theta) H Rz(lam - pi/2).
    """
    theta, phi, lam = theta_phi_lambda
    
    def synthesize_rz(angle, q):
        if abs(angle) < 1e-15:
            return []
        circuit = pygridsynth.gridsynth_circuit(mp.mpf(angle), mp.mpf("1e-12"))
        qasm = []
        for gate in circuit:
            name = type(gate).__name__
            if name == "HGate": qasm.append(f"h q[{q}];")
            elif name == "TGate": qasm.append(f"t q[{q}];")
            elif name == "SGate": qasm.append(f"t q[{q}]; t q[{q}];")
            elif name == "SXGate":
                qasm.append(f"h q[{q}]; t q[{q}]; t q[{q}]; t q[{q}]; t q[{q}]; h q[{q}];")
        return qasm

    # Decomposition: U3(theta, phi, lam) = Rz(phi + pi/2) H Rz(theta) H Rz(lam - pi/2)
    gates = synthesize_rz(lam - np.pi/2, qubit)
    gates.append(f"h q[{qubit}];")
    gates.extend(synthesize_rz(theta, qubit))
    gates.append(f"h q[{qubit}];")
    gates.extend(synthesize_rz(phi + np.pi/2, qubit))
    return gates

def synthesize_unitary_2q(U, filename):
    """
    Transpiles a 2-qubit unitary and synthesizes it into Clifford+T gates.
    """
    qc = QuantumCircuit(2)
    qc.append(Operator(U), [0, 1])
    qc_decomposed = transpile(qc, basis_gates=["cx", "u3"], optimization_level=3)
    
    qasm_lines = ["OPENQASM 2.0;", "include \"qelib1.inc\";", "qreg q[2];"]
    
    for instr, qargs, cargs in qc_decomposed.data:
        if instr.name == "cx":
            c = qc_decomposed.find_bit(qargs[0]).index
            t = qc_decomposed.find_bit(qargs[1]).index
            qasm_lines.append(f"cx q[{c}], q[{t}];")
        elif instr.name in ["u3", "u"]:
            params = [float(p) for p in instr.params]
            q = qc_decomposed.find_bit(qargs[0]).index
            qasm_lines.extend(get_clifford_t_qasm(params, q))
            
    with open(filename, "w") as f:
        f.write("\n".join(qasm_lines) + "\n")

def get_h3_unitary(t):
    # Task 6 Hamiltonian H3 = XX + ZI + IZ
    X = np.array([[0, 1], [1, 0]])
    Z = np.array([[1, 0], [0, -1]])
    I = np.eye(2)
    XX = np.kron(X, X)
    ZI = np.kron(Z, I)
    IZ = np.kron(I, Z)
    H = XX + ZI + IZ
    return la.expm(1j * t * H)

def solve_all():
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    
    # Task 3: exp(i pi/7 ZZ)
    theta3 = np.pi/7
    U3 = la.expm(1j * theta3 * np.kron(Z, Z))
    synthesize_unitary_2q(U3, "task3_expzz.qasm")
    
    # Task 4: exp(i pi/7 (XX + YY))
    theta4 = np.pi/7
    U4 = la.expm(1j * theta4 * (np.kron(X, X) + np.kron(Y, Y)))
    synthesize_unitary_2q(U4, "task4_h1.qasm")
    
    # Task 6: exp(i pi/7 (XX + ZI + IZ))
    theta6 = np.pi/7
    U6 = get_h3_unitary(theta6)
    synthesize_unitary_2q(U6, "task6_h3.qasm")
    
    # Task 7: State preparation
    import qiskit.quantum_info as qi
    target_state = qi.random_statevector(4, seed=42).data
    from qiskit.circuit.library import StatePreparation
    sp = StatePreparation(target_state)
    qc7 = transpile(sp.definition, basis_gates=["cx", "u3"], optimization_level=3)
    synthesize_unitary_2q(Operator(qc7).data, "task7_state.qasm")
    
    # Task 8: Structured 1
    U8 = 0.5 * np.array([
        [1, 1, 1, 1],
        [1, 1j, -1, -1j],
        [1, -1, 1, -1],
        [1, -1j, -1, 1j]
    ])
    synthesize_unitary_2q(U8, "task8_u1.qasm")
    
    # Task 9: Structured 2
    a = (1 - 1j) / 2
    b = (1 + 1j) / 2
    U9 = np.array([
        [1, 0, 0, 0],
        [0, 0, -a, b],
        [0, 1j, 0, 0],
        [0, 0, -a, -b]
    ])
    synthesize_unitary_2q(U9, "task9_u2.qasm")
    
    # Task 10: Random Unitary
    U10 = qi.random_unitary(4, seed=42).data
    synthesize_unitary_2q(U10, "task10_random.qasm")

    # Task 11: 4-qubit Diagonal
    phases_idx = [0, 4, 6, 5, 7, 7, 6, 5, 5, 7, 6, 6, 6, 6, 7, 5]
    actual_phases = [p * np.pi / 4 for p in phases_idx]
    amplitudes = [np.exp(1j * p) for p in actual_phases]
    from qiskit.circuit.library import Diagonal
    qc11 = Diagonal(amplitudes)
    qc11_decomposed = transpile(qc11, basis_gates=["cx", "u3", "t", "tdg", "s", "sdg", "z"], optimization_level=3)
    
    qasm11 = ["OPENQASM 2.0;", "include \"qelib1.inc\";", "qreg q[4];"]
    for instr, qargs, cargs in qc11_decomposed.data:
        q_indices = [qc11_decomposed.find_bit(q).index for q in qargs]
        if instr.name == "cx":
            qasm11.append(f"cx q[{q_indices[0]}], q[{q_indices[1]}];")
        elif instr.name == "t":
            qasm11.append(f"t q[{q_indices[0]}];")
        elif instr.name == "tdg":
            qasm11.append(f"tdg q[{q_indices[0]}];")
        elif instr.name == "s":
            qasm11.append(f"t q[{q_indices[0]}]; t q[{q_indices[0]}];")
        elif instr.name == "sdg":
            qasm11.append(f"tdg q[{q_indices[0]}]; tdg q[{q_indices[0]}];")
        elif instr.name == "z":
            qasm11.append(f"t q[{q_indices[0]}]; t q[{q_indices[0]}]; t q[{q_indices[0]}]; t q[{q_indices[0]}];")
        elif instr.name == "u3":
            p = [float(x) for x in instr.params]
            qasm11.extend(get_clifford_t_qasm(p, q_indices[0]))
            
    with open("task11_diagonal.qasm", "w") as f:
        f.write("\n".join(qasm11) + "\n")

if __name__ == "__main__":
    solve_all()
