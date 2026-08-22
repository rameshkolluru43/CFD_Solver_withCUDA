import os 
import numpy as np
from numpy import complex128
import numpy.typing as npt
import matplotlib as plt
import scipy.linalg as spla
from scipy.linalg import norm, inv
from scipy.optimize import minimize
from scipy.spatial import distance
import scipy
import math
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import Operator
#from qiskit_aer import Aer, QasmSimulator
from qiskit_algorithms.optimizers import COBYLA, SPSA
from qiskit.compiler import transpile
from qiskit.circuit.library import TwoLocal
import matplotlib.pyplot as plt
from typing import List, TypeVar, Union, cast 


ATOL_DEFAULT = 1e-8
RTOL_DEFAULT = 1e-5
complex_type = TypeVar("complex_type", float, complex)
complex_array_type = npt.NDArray[np.cdouble]

# Using this Exception class to avoid more iterations than needed, and to stop at the point that we define
class MaxIterationsReached(Exception):
    pass

class Utility: 
    def __init__(self):
        cwd = os.getcwd() 
        print(f"Current working directory: {cwd}\n")


    def normalize(self, v: complex_array_type) -> complex_array_type:
        """Normalize the given vector."""
        norm_arr = scipy.linalg.norm(v)
        return v / norm_arr, norm_arr 

    def normalize_matrix(self, A: complex_array_type) -> complex_array_type:
        """Normalize the given matrix."""
        det = np.linalg.det(A)
        return A / det, det

    def normL2norm(self, qu, cl): 
        """Calculation of L2 Norm between quantum generated vector and classical vector"""
        l2_norm = np.sum(np.power((qu-cl),2)) #previously normalize_array(cl)
        y = np.sum(np.power((cl),2)) 
        res = np.sqrt(l2_norm/y)
        return res


    def convert_bin(self, num, ancilla_size= 1): 
        '''Function changes a decimal number to binary number of deired length. By default the length is 1'''
        n = ancilla_size 
        b = bin(num)[2:]
        l = len(b)
        b = str(0) * (n - l) + b
        b = b[::-1]
        return b

    def auxilliary_matrix(self, 
            x: Union[npt.NDArray[np.float64], complex_array_type]
        ) -> complex_array_type:
            """Returns the auxiliary matrix for the decomposition of size n.
            derived and defined as : i * sqrt(I - x^2)

            Args:
                x (Union[npt.NDArray[np.float_], complex_array_type]): original matrix.

            Returns:
                complex_array_type: The auxiliary matrix.
            """
            mat = np.eye(len(x)) - x @ x
            mat = cast(npt.NDArray[Union[np.float64, np.cdouble]], spla.sqrtm(mat))
            return 1.0j * mat
    
    def decompose(self, A: complex_array_type, b: complex_array_type) -> complex_array_type:
        """Decompose a given matrix A into B and C.

        Args:
            A (complex_array_type): The matrix to be decomposed.
            b (complex_array_type): The vector to be normalized. 

        Returns:
            Ub, Vb, Uc, Vc: The decomposed matrices, it's coefficients, and normlized b vector. 
        """
        global A_normalized, b_normalized, norm_A, norm_b
        A, norm_A = self.normalize(A)
        A_normalized = A 
        b, norm_b = self.normalize(b)
        b_normalized = b
        B = 0.5*(A+A.T)
        C = (0.5/1j)*(A-A.T)

        # hold the Ub,Vb,Uc,Vc matrices in unitary_matrices and coefficents in the unitary_coefficients.
        unitary_matrices, unitary_coefficients = [], [] 

        # calcualation of coefficents.
        coef_real = norm_A*0.5
        coef_imag = coef_real * 1j

        # check for norm_B and norm_C <=1
        print("cheking for norm_B and norm_C <=1")
        norm_B = np.linalg.norm(B) <= 1
        norm_C = np.linalg.norm(C) <= 1
        print(f"Norm_B: {norm_B}, Norm_C: {norm_C}\n")

        #creation of Ub,Vb unitary matrices
        if not np.allclose(B, 0.0):
                    aux_mat = self.auxilliary_matrix(B)
                    unitary_matrices += [B + aux_mat, B - aux_mat] 
                    unitary_coefficients += [coef_real] * 2

        # creation Uc,Vc matrices
        if not np.allclose(C, 0.0):
                    aux_mat = self.auxilliary_matrix(C)
                    unitary_matrices += [C + aux_mat, C - aux_mat]
                    unitary_coefficients += [coef_imag] * 2

        unit_coeffs = np.array(unitary_coefficients, dtype=np.cdouble)
        return unitary_matrices, unit_coeffs, norm_b
    
    def is_identity_matrix(self, mat, ignore_phase=False, rtol=RTOL_DEFAULT, atol=ATOL_DEFAULT):
        """Test if an array is an identity matrix."""
        if atol is None:
            atol = ATOL_DEFAULT
        if rtol is None:
            rtol = RTOL_DEFAULT
        mat = np.array(mat)
        if mat.ndim != 2:
            return False
        if ignore_phase:
            """If the matrix is equal to an identity up to a phase, we can
            remove the phase by multiplying each entry by the complex
            conjugate of the phase of the [0, 0] entry. """
            theta = np.angle(mat[0, 0])
            mat = np.exp(-1j * theta) * mat
        # Check if square identity
        iden = np.eye(len(mat))
        return np.allclose(mat, iden, rtol=rtol, atol=atol)

    def is_A_matrix_accurate(self, unitaries: List) -> float:
        """ to evaluate the given matrices representation accuracy.

        Args:
            unitaries (List): 4 unitary matrix list.

        Returns:
            float: returns a the Forbenius norm of the difference between the original and reconstructed matrix.
        """
        if len(unitaries) != 4:
            A_reconstructed = 0.5*(unitaries[0]+unitaries[1])
        else:
            A_reconstructed = 0.5*(unitaries[0]+unitaries[1]+unitaries[2]+unitaries[3])
        error = np.linalg.norm(A_normalized - A_reconstructed, 'fro')
        trace_A_reconstructed = np.trace(A_reconstructed)
        trace_A_normalized = np.trace(A_normalized)
        print(f"Trace of A_reconstructed: {trace_A_reconstructed}")
        print(f"Trace of A_normalized: {trace_A_normalized}")
        print(f"Reconstruction error (Frobenius norm): {error}\n")
        
    def define_A_b(self, A_pass, b_pass): 
        if not isinstance(A_pass, np.ndarray):
            A = np.array(A_pass)
        else:
            A = A_pass
        if not isinstance(b_pass, np.ndarray):
            b = np.array(b_pass)
        else:
            b = b_pass
        
        if b.ndim != 1:
            raise ValueError("b_pass must be a 1-dimensional vector.")
        si = int(A[0].size)
        Ni = math.ceil(math.log(len(A[0]),2))
        N=int(math.pow(2,math.ceil(math.log(A[0].size,2))))-si
        N_cut = N
        qq = np.zeros((2**(Ni), 2**(Ni)))
        for i in range(si):
            for j in range(si):
                qq[i,j] = A[i,j]
        # for k in range(Ni): 
        #     for 
        
        eps = 0.00001
        Iq = np.eye(2**Ni,2**Ni)
        A = qq + (Iq*eps)
        b = np.pad(b, ((0,N)), mode='constant')
        print(A.shape)
        return A, b, N_cut


class Qsolver(Utility,MaxIterationsReached): 
    def __init__(self, A, b):
        # self.A = A 
        # self.b = b
        self.A, self.b, self.N_cut = self.define_A_b(A,b)
        # Decomposing the Unitaary matrices, and checking the accuracy of the reconstructed matrix.
        self.unitary_matrices, self.unitary_coefficients, self.norm_b = self.decompose(self.A, self.b)
        print("Decomposition Completed")
        self.ancilla_size = 2 
        self.system_size = int(np.log2(self.A.shape[0])) 
        self.tot_size = self.system_size + self.ancilla_size 
        self.iteration_count = 0
        self.max_iterations = 2
        self.opt_param = None
        self.reps = 6
        self.cost_values = []
        self.w_init = [np.pi/2] * (self.reps * ((self.system_size*2)-2)+self.system_size)
        self.opt_param = None
        # self.detu = scipy.linalg.det(self.A)
        # print("Evaluating the accuracy of the reconstructed matrix ")
        # self.is_A_matrix_accurate(self.unitary_matrices)
        # print(f"Checking the determinant of the matrix A: {self.detu}")
        # self.mat_det = scipy.linalg.det(self.A/self.detu)
        # print(f"Checking the determinant of the normlized matrix A: {self.mat_det}")
        # print()
        print("non-singular matrix check: ")
        print(self.A.dot(inv(self.A)))
        print("Condition number of the matrix A: ") 
        print(np.linalg.cond(self.A))
        # checking for U*(U^H) = I and V*(V^H) = I
        print("checking for U*(U^H) = I and V*(V^H) = I\n")
        if len(self.unitary_matrices) != 4:
            Ub_is_unitary = np.conj(self.unitary_matrices[0].T).dot(self.unitary_matrices[0])
            Vb_is_unitary = np.conj(self.unitary_matrices[1].T).dot(self.unitary_matrices[1])
            print(self.is_identity_matrix(Ub_is_unitary))
            print(self.is_identity_matrix(Vb_is_unitary))
        else:
            Ub_is_unitary = np.conj(self.unitary_matrices[0].T).dot(self.unitary_matrices[0])
            Vb_is_unitary = np.conj(self.unitary_matrices[1].T).dot(self.unitary_matrices[1])
            Uc_is_unitary = np.conj(self.unitary_matrices[2].T).dot(self.unitary_matrices[2])
            Vc_is_unitary = np.conj(self.unitary_matrices[3].T).dot(self.unitary_matrices[3])

            print(self.is_identity_matrix(Ub_is_unitary))
            print(self.is_identity_matrix(Vb_is_unitary))
            print(self.is_identity_matrix(Uc_is_unitary))
            print(self.is_identity_matrix(Vc_is_unitary))


    def create_unitary_circuits(self,
            unimatrices: List[np.ndarray], names: List[str]
        ) -> List[QuantumCircuit]:
            """Construct the quantum circuits from unitary matrices

            Args:
                unimatrices (List[np.ndarray]): list of unitary matrices of the decomposition.
                names (List[str]): names of the circuits

            Returns:
                List[QuantumCircuit]: quantum circuits
            """

            def make_qc(mat: complex_array_type, name: str) -> QuantumCircuit:
                circuit = QuantumCircuit(self.system_size, name=name)
                circuit.unitary(mat, circuit.qubits)
                return circuit

            return [make_qc(mat, name) for mat, name in zip(unimatrices, names)]


    def add_controls(self, num_controls: int) -> QuantumCircuit:
        """Add control qubits to the unitary circuit.

        Args:
            num_controls (int): The number of control qubits to add.

        Returns:
            QuantumCircuit: The circuit with added control qubits.
        """
        Q_circuit = QuantumCircuit(self.tot_size)
        Q_circuit.h(range(self.ancilla_size)) 
        names = ["Ub", "Vb","Uc","Vc"]
        circuits = self.create_unitary_circuits(self.unitary_matrices, names)
        if len(circuits) != 4:
            controlled_U1 = circuits[0].control(num_controls, ctrl_state="00")
            controlled_U2 = circuits[1].control(num_controls, ctrl_state="01")
            controlled_circ_lst = [controlled_U1, controlled_U2]
            Q_circuit.append(controlled_circ_lst[0].to_instruction(), range(self.tot_size))
            Q_circuit.append(controlled_circ_lst[1].to_instruction(), range(self.tot_size))
            return Q_circuit
        else:
            controlled_U1 = circuits[0].control(num_controls, ctrl_state="00")  
            controlled_U2 = circuits[1].control(num_controls, ctrl_state="01")
            controlled_U3 = circuits[2].control(num_controls, ctrl_state="10")
            controlled_U4 = circuits[3].control(num_controls, ctrl_state="11")
            controlled_circ_lst = [controlled_U1, controlled_U2, controlled_U3, controlled_U4]

            Q_circuit.append(controlled_circ_lst[0].to_instruction(), range(self.tot_size))
            Q_circuit.append(controlled_circ_lst[1].to_instruction(), range(self.tot_size))
            Q_circuit.append(controlled_circ_lst[2].to_instruction(), range(self.tot_size))
            Q_circuit.append(controlled_circ_lst[3].to_instruction(), range(self.tot_size))

            return Q_circuit
            
    def variational_ansatz(self, parameters: List[float]) -> QuantumCircuit:
        """Construct the variational ansatz circuit.
        
        Args:
            parameters (List[float]): The parameters for the circuit.

        Returns:
            QuantumCircuit: The variational ansatz circuit.
        """
        global para_count

        circuit = QuantumCircuit(self.system_size)

        # alternating Controlled-Z operations
        if self.system_size > 3:
            if self.system_size % 2 == 0:
                for i in range(0, self.system_size, 2):
                    circuit.cz(i, i + 1)
            else: 
                for i in range(0, self.system_size-1, 2):
                    circuit.cz(i, i + 1)
        else:
            circuit.cz(0, 1)


        # Ry gates to all qubits
        para_count = 0
        for i in range(self.system_size): 
            circuit.ry(parameters[i], i) 
            para_count += 1

        # Controlled-Z operations alternating, but not for the first and last qubits.
        if self.system_size > 3:
            for i in range(1, self.system_size-1, 2):
                circuit.cz(i, i + 1)
        else:
            circuit.cz(1, 2)

        # Ry gates except the first and last qubits.
        for i in range(1, self.system_size-1):
            circuit.ry(parameters[para_count], i)
            para_count += 1
        return circuit 

    def variational_circuit(self, reps: int, parameters: List[float]) -> QuantumCircuit:
        """Construct the variational circuit.

        Args:
            reps (int): The number of repetitions.
            parameters (List[float]): The parameters for the circuit.

        Returns:
            QuantumCircuit: The variational circuit.
        """

        circuit = QuantumCircuit(self.system_size, name="V(w)")

        qbit_per_ansatz = (self.system_size*2)-2
        for i in range(self.system_size):
            circuit.ry(parameters[i], i)
        parameters = parameters[self.system_size:]
        for _ in range(reps):
            circuit.append(self.variational_ansatz(parameters[:qbit_per_ansatz]), range(self.system_size))
            parameters = parameters[qbit_per_ansatz:]

        return circuit

    def embedding_b(self, b) -> QuantumCircuit:
        """Embed the vector b into the quantum circuit.

        Args:
            b : The vector to embed.

        Returns:
            QuantumCircuit: The embedded vector.
        """

        b_vector_circuit = QuantumCircuit(self.system_size,name="bvec")
        b_vector_circuit.prepare_state(b)
        return b_vector_circuit


    def full_circuit(self, parameters: List[float]) -> QuantumCircuit:
        """Construct the full quantum circuit.""" 

        print(self.tot_size)
        F_circuit = QuantumCircuit(self.tot_size) 
        F_circuit.append(self.add_controls(self.ancilla_size), range(self.tot_size)) 
        F_circuit.append(self.variational_circuit(self.reps, parameters), range(2,self.tot_size))
        F_circuit.append(self.embedding_b(b_normalized), range(2,self.tot_size)) 
        F_circuit.measure_all()
        return F_circuit

    def prepare_x(self, opt_weights: List[float], reps: int) -> np.ndarray:
        """Prepare the x vector"""
        print("entered prep_x fun")
        Var_circuit = QuantumCircuit(self.system_size)
        Var_circuit.h(range(self.system_size))
        Var_circuit.append(self.variational_circuit(reps, opt_weights), range(self.system_size))  


        # Sampling 
        Var_circuit.measure_all()
        backend = QasmSimulator() 
        qc_prepx = transpile(Var_circuit, backend=backend)
        job = backend.run(qc_prepx, shots=1024)
        result = job.result()
        counts = result.get_counts()
        old_keys = [self.convert_bin(_, self.system_size) for _ in range(2**self.system_size)]
        new_keys = [old_keys[i][::-1] for i in range(2**self.system_size)]
        new_count = {new_keys[i]: counts.get(new_keys[-(i+1)], 0) for i in range(2**self.system_size)}

        stv = np.array(list(new_count.values()))
        x = stv / 10**4 
        print("completed preparing x")
        return x
        


    def optimize_cost_function(self, optimizer: str):
        global iteration_count
        iteration_count = 0
        
        options = {
            'maxiter': self.max_iterations,
            'disp': True,
            'rhobeg': 1.0,
            'tol': 1e-7
        }
        if optimizer == 'COBYLA':
            optimizer = COBYLA(options=options)
        elif optimizer == 'SPSA':
            print("Using SPSA optimizer")
            optimizer = SPSA(maxiter=self.max_iterations,learning_rate=0.5,perturbation=0.999)
        
        try:
            result = optimizer.minimize(self.cost, self.w_init,bounds=[(0, 2 * np.pi) for i in range(len(self.w_init))])
            self.opt_param = result.x 
        except MaxIterationsReached:
            print(f"Maximum iterations of {self.max_iterations} reached.")
        
        return self.opt_param

    def cost(self, weights):
        print("entered cost fn")
        print(self.iteration_count)
        self.iteration_count += 1
        if self.iteration_count > 10:
            raise MaxIterationsReached
        
        # Will be using it if we want to store the best parameters among all iterations
        # if opt_param is None or cost_value < cost(opt_param ):
        #   opt_param = weights

        self.opt_param  = weights
        print(weights)
        print(self.opt_param)
        qc = self.full_circuit(self.opt_param)
        print(qc)
        backend = QasmSimulator() 
        print(backend)
        qc_cost = transpile(qc, backend=backend)
        print(qc_cost)
        job = backend.run(qc_cost, shots=1024)
        print(job)
        result = job.result()
        print(result)
        counts = result.get_counts()
        print(counts)

        old_keys = [self.convert_bin(_, self.tot_size) for _ in range(2**self.tot_size)]
        new_keys = [old_keys[i][::-1] for i in range(2**self.tot_size)]
        new_count = {new_keys[i]: counts.get(new_keys[-(i+1)], 0) for i in range(2**self.tot_size)}

        stv = np.array(list(new_count.values()))
        statevector = stv / 10**4 

        P_global = np.zeros((2**self.tot_size, 2**self.tot_size))
        P_global[0, 0] = 1.0

        P_ancilla = np.zeros((2**self.ancilla_size, 2**self.ancilla_size)) 
        P_ancilla[0, 0] = 1.0

        P_ancilla_full = np.kron(np.eye(2**self.system_size), P_ancilla)
        proj_global = Operator(P_global)
        proj_ancilla = Operator(P_ancilla_full)

        expval_global = np.real(np.dot(statevector.conj(), np.dot(proj_global, statevector)))
        expval_ancilla = np.real(np.dot(statevector.conj(), np.dot(proj_ancilla, statevector)))

        p_cond = expval_global / expval_ancilla
        cost_value = 1 - p_cond
        self.cost_values.append(cost_value)
        print(f"Cost: {cost_value}")
        return cost_value

    def fine_tune_cost(self, weight , x):
        '''
        Cost function for Fine tuning.
        '''
        weighted_x = np.array([weight[i]*x[i] for i in range(len(self.b))])
        b_made = np.matmul(self.A,weighted_x)
        cost = distance.euclidean(b_made ,self.b)
        return(cost)

    def Fine_Tune(self, x):
        '''
        [x1,x2,x3...] ---> [k1x1,k2x2,k3x3...]
        k1 ,k2 , k3 are estimated using an optimiser.
        '''
        weights = np.ones(len(self.b))
        weights = minimize(self.fine_tune_cost,weights,args=(x))
        sol = np.array([weights['x'][i]*x[i] for i in range(len(x))])
        return sol

    def solve(self, optimizer: str,callFnT: bool) -> np.ndarray: 
        "This is the main function"
        print("entered solve")
        opt_weights = self.optimize_cost_function(optimizer)
        x_vector = self.prepare_x(opt_weights,self.reps)
        print("comleted preperation of x")
        if callFnT == True: 
            print("Entered Fine Tune")
            fin_sol = self.Fine_Tune(x_vector)
            print(fin_sol[:self.N_cut])
            return fin_sol[:self.N_cut]
        else:
            print(x_vector[:self.N_cut])
            return x_vector[:self.N_cut]
        
