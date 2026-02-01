import numpy as np
import qiskit
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, state_fidelity
import sys

def get_distance(U_target, U_compiled):
    """
    Calculates the operator norm distance:
    d(U, U_tilde) = min_{phi \in [0, 2pi)} ||U - e^{i phi} U_tilde||_op
    """
    # The operator norm is the largest singular value of (U - e^{i phi} U_tilde)
    # This is equivalent to min_phi || I - e^{i phi} U^\dagger U_tilde ||_op
    
    # Let A = U^\dagger @ U_compiled
    A = np.conj(U_target).T @ U_compiled
    
    # We want to find phi that minimizes || U - e^{i phi} U_tilde ||
    # This is equivalent to finding phi that makes e^{i phi} * eigvals(A) as 
    # close to 1 as possible. 
    # Actually, the challenge says 'operator norm distance'.
    # A common way to handle global phase is to align the phase of the trace
    # or to look at the eigenvalues.
    
    # For unitary matrices, the operator norm of (U - e^{i phi} U_tilde) 
    # is related to the eigenvalues of U^\dagger U_tilde.
    # Let the eigenvalues of A be exp(i theta_j).
    # then || U - e^{i phi} U_tilde || = max_j | 1 - exp(i (theta_j - phi)) |
    
    eigvals = np.linalg.eigvals(A)
    phases = np.angle(eigvals)
    
    # To minimize the max distance, we align with the first eigenvalue
    # and then possibly center the remaining spread.
    shifted_phases = (phases - phases[0] + np.pi) % (2 * np.pi) - np.pi
    
    # Optional: centering for even better precision
    mid_shift = (np.min(shifted_phases) + np.max(shifted_phases)) / 2
    shifted_phases -= mid_shift
    
    dist = np.max(np.abs(1 - np.exp(1j * shifted_phases)))
    return dist

def get_t_count(circuit):
    t_count = 0
    for instruction in circuit.data:
        name = instruction.operation.name
        if name in ['t', 'tdg']:
            t_count += 1
    return t_count

def get_qiskit_unitary(target_2x2, control=0, target=1, num_qubits=2):
    """
    Constructs a controlled unitary in Qiskit ordering.
    In Qiskit, for a 2-qubit system:
    Basis is |q1 q0>: |00>, |01>, |10>, |11>
    If q0 is control and q1 is target:
    |00> -> |00>
    |01> -> |1 q0> -> |1 1> if target is applied? No.
    Let's use Operator.from_label or similar if possible, 
    but for now, explicit construction.
    """
    from qiskit.quantum_info import Operator
    from qiskit import QuantumCircuit
    qc = QuantumCircuit(num_qubits)
    # This is a bit circular, but we can use qiskit to build the 'exact' one
    # using high-level gates, then compare with Clifford+T.
    return Operator(qc) # Placeholder

def verify_unitary(qasm_file, target_unitary, verbose=True):
    """
    target_unitary should be a numpy array in Qiskit ordering.
    """
    qc = QuantumCircuit.from_qasm_file(qasm_file)
    compiled_unitary = Operator(qc).data
    
    if isinstance(target_unitary, Operator):
        target_unitary = target_unitary.data

    dist = get_distance(target_unitary, compiled_unitary)
    t_count = get_t_count(qc)
    
    if verbose:
        print(f"File: {qasm_file}")
        print(f"Operator Norm Distance: {dist:.10e}")
        print(f"T-count: {t_count}")
        
        allowed = {'h', 't', 'tdg', 'cx'}
        for instruction in qc.data:
            if instruction.operation.name not in allowed:
                print(f"WARNING: Forbidden gate '{instruction.operation.name}' used!")
            
    return dist, t_count

if __name__ == "__main__":
    # Example usage:
    # target = np.array([[1, 0], [0, 1]])
    # verify_unitary("test.qasm", target)
    pass
