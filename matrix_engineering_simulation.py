#!/usr/bin/env python3
"""
===============================================================================
ENGINEERING MATRIX SOLUTIONS SIMULATION & TASK SUITE
===============================================================================
A task-based Python engineering matrix simulation designed for students and 
executable on GitHub (CLI, GitHub Actions CI/CD Auto-grader, Interactive mode).

Includes three real-world engineering matrix systems:
  Task 1: The Mechanical Truss (Structural Engineering - Joint Equilibrium Matrix)
  Task 2: The Electrical Circuit (Electrical Engineering - Nodal Analysis Matrix)
  Task 3: The Fluid Mixing Tanks (Chemical Engineering - Steady-State Mass Balance Matrix)

Author: Gemini Notebook Educational Series
Target Audience: HNC / Undergraduate Engineering Students (Level 4/5)
===============================================================================
"""

import sys
import math
import json
import argparse
from typing import List, Tuple, Dict, Any

# -----------------------------------------------------------------------------
# PURE PYTHON MATRIX UTILITIES (Zero External Dependency Fallback)
# -----------------------------------------------------------------------------
class MatrixOps:
    """Pure Python matrix linear algebra routines for educational transparent execution."""
    
    @staticmethod
    def create_zeros(rows: int, cols: int) -> List[List[float]]:
        return [[0.0 for _ in range(cols)] for _ in range(rows)]

    @staticmethod
    def mat_mul(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
        r_a, c_a = len(A), len(A[0])
        r_b, c_b = len(B), len(B[0])
        if c_a != r_b:
            raise ValueError(f"Cannot multiply matrix {r_a}x{c_a} with {r_b}x{c_b}")
        C = MatrixOps.create_zeros(r_a, c_b)
        for i in range(r_a):
            for j in range(c_b):
                C[i][j] = sum(A[i][k] * B[k][j] for k in range(c_a))
        return C

    @staticmethod
    def mat_vec_mul(A: List[List[float]], x: List[float]) -> List[float]:
        r, c = len(A), len(A[0])
        if c != len(x):
            raise ValueError(f"Matrix cols ({c}) != Vector length ({len(x)})")
        return [sum(A[i][j] * x[j] for j in range(c)) for i in range(r)]

    @staticmethod
    def determinant(A: List[List[float]]) -> float:
        n = len(A)
        if n != len(A[0]):
            raise ValueError("Determinant requires a square matrix")
        if n == 1:
            return A[0][0]
        if n == 2:
            return A[0][0] * A[1][1] - A[0][1] * A[1][0]
        
        # Gaussian Elimination for Det calculation
        M = [row[:] for row in A]
        det = 1.0
        for i in range(n):
            pivot = i
            for j in range(i + 1, n):
                if abs(M[j][i]) > abs(M[pivot][i]):
                    pivot = j
            if pivot != i:
                M[i], M[pivot] = M[pivot], M[i]
                det *= -1.0
            if abs(M[i][i]) < 1e-12:
                return 0.0
            det *= M[i][i]
            for j in range(i + 1, n):
                factor = M[j][i] / M[i][i]
                for k in range(i, n):
                    M[j][k] -= factor * M[i][k]
        return det

    @staticmethod
    def solve_gaussian(A: List[List[float]], b: List[float]) -> List[float]:
        """Solves A x = b using Gaussian Elimination with partial pivoting."""
        n = len(A)
        # Augment matrix A with vector b
        M = [A[i][:] + [b[i]] for i in range(n)]
        
        # Forward elimination
        for i in range(n):
            pivot = i
            for j in range(i + 1, n):
                if abs(M[j][i]) > abs(M[pivot][i]):
                    pivot = j
            if pivot != i:
                M[i], M[pivot] = M[pivot], M[i]
            
            if abs(M[i][i]) < 1e-12:
                raise ValueError("Matrix is singular or near-singular. Cannot solve A x = b.")
            
            for j in range(i + 1, n):
                factor = M[j][i] / M[i][i]
                for k in range(i, n + 1):
                    M[j][k] -= factor * M[i][k]
                    
        # Back substitution
        x = [0.0] * n
        for i in range(n - 1, -1, -1):
            s = sum(M[i][j] * x[j] for j in range(i + 1, n))
            x[i] = (M[i][n] - s) / M[i][i]
        return x

    @staticmethod
    def invert(A: List[List[float]]) -> List[List[float]]:
        """Inverts square matrix A using Gauss-Jordan Elimination."""
        n = len(A)
        # Identity matrix augmented
        M = [A[i][:] + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        
        for i in range(n):
            pivot = i
            for j in range(i + 1, n):
                if abs(M[j][i]) > abs(M[pivot][i]):
                    pivot = j
            if pivot != i:
                M[i], M[pivot] = M[pivot], M[i]
            
            p_val = M[i][i]
            if abs(p_val) < 1e-12:
                raise ValueError("Matrix is non-invertible (singular).")
            
            # Normalize row i
            for k in range(2 * n):
                M[i][k] /= p_val
                
            # Eliminate column i in other rows
            for j in range(n):
                if j != i:
                    factor = M[j][i]
                    for k in range(2 * n):
                        M[j][k] -= factor * M[i][k]
                        
        inv = [[M[i][j + n] for j in range(n)] for i in range(n)]
        return inv


# -----------------------------------------------------------------------------
# TASK 1: THE MECHANICAL TRUSS SIMULATION
# -----------------------------------------------------------------------------
class MechanicalTrussTask:
    """
    Structural Engineering Task:
    A pin-jointed 2D bridge truss subjected to downward force P = 60 kN at Joint 2.
    """
    
    def __init__(self):
        self.title = "Task 1: The Mechanical Truss Equilibrium"
        self.P = 60.0  # kN applied vertical load at Joint 2
        self.angle = 45.0 # degrees
        self.rad = math.radians(self.angle)
        self.cos45 = round(math.cos(self.rad), 4)
        self.sin45 = round(math.sin(self.rad), 4)
        
        self.A = [
            [ self.cos45,  0.0,         1.0,  1.0,  0.0],
            [ self.sin45,  0.0,         0.0,  0.0,  1.0],
            [-self.cos45,  self.cos45,  0.0,  0.0,  0.0],
            [-self.sin45, -self.sin45,  0.0,  0.0,  0.0],
            [ 0.0,        -self.cos45, -1.0,  0.0,  0.0]
        ]
        self.b = [0.0, 0.0, 0.0, -self.P, 0.0]
        self.variable_names = ["f1 (Member 1)", "f2 (Member 2)", "f3 (Member 3)", "R1x (Rx pin)", "R1y (Ry pin)"]

    def display_ascii(self):
        return r"""
       MECHANICAL TRUSS SCHEMATIC (Joint Forces):
                 (Joint 2)
                   / \
           m1     /   \     m2
                 /  P  \  (P = 60 kN down)
                /   v   \
      (Joint 1)-----------(Joint 3)
         /\       m3        O (Roller)
      (Pin R1x, R1y)     (R3y)
        """

    def solve(self) -> List[float]:
        return MatrixOps.solve_gaussian(self.A, self.b)


# -----------------------------------------------------------------------------
# TASK 2: THE ELECTRICAL CIRCUIT SIMULATION
# -----------------------------------------------------------------------------
class ElectricalCircuitTask:
    """
    Electrical Engineering Task:
    Nodal Analysis of a 3-Node DC Bridge Circuit with resistors and power sources.
    """

    def __init__(self):
        self.title = "Task 2: Electrical Nodal Conductance Matrix"
        self.G = [
            [ 0.75, -0.25,  0.0 ],
            [-0.25,  0.95, -0.20],
            [ 0.00, -0.20,  0.30]
        ]
        self.I = [12.0, 0.0, 4.0]
        self.variable_names = ["V1 (Volts)", "V2 (Volts)", "V3 (Volts)"]

    def display_ascii(self):
        return r"""
       ELECTRICAL CIRCUIT NODAL SCHEMATIC:
       (+12A) ---> (V1) ----[R2=4Ω]---- (V2) ----[R4=5Ω]---- (V3) <--- (+4A)
                    |                    |                    |
                 [R1=2Ω]              [R3=2Ω]              [R5=10Ω]
                    |                    |                    |
                  (GND)                (GND)                (GND)
        """

    def solve(self) -> List[float]:
        return MatrixOps.solve_gaussian(self.G, self.I)

    def get_impedance_matrix(self) -> List[List[float]]:
        """Returns the inverse matrix Z = G^(-1) representing the nodal resistance matrix."""
        return MatrixOps.invert(self.G)


# -----------------------------------------------------------------------------
# TASK 3: THE FLUID MIXING TANKS SIMULATION
# -----------------------------------------------------------------------------
class FluidMixingTanksTask:
    """
    Chemical & Process Engineering Task:
    Continuous Stirred-Tank Reactor (CSTR) system with 3 interconnected fluid mixing tanks.
    """

    def __init__(self):
        self.title = "Task 3: Fluid Mixing Tanks Steady-State Mass Balance"
        self.M = [
            [ 10.0,   0.0,  -3.0],
            [-10.0,  10.0,   0.0],
            [  0.0,  -8.0,  13.0]
        ]
        self.b = [500.0, 0.0, 100.0]
        self.variable_names = ["c1 (Tank 1 g/L)", "c2 (Tank 2 g/L)", "c3 (Tank 3 g/L)"]

    def display_ascii(self):
        return r"""
       FLUID MIXING TANKS SCHEMATIC:
       Inlet 1 (10 L/min, 50 g/L)
            |
            v
       +----------+    Q12 = 10 L/min   +----------+
       |  TANK 1  | ------------------> |  TANK 2  | ----> Exit 2 (2 L/min)
       +----------+                     +----------+
            ^                                |
            | Q31 = 3 L/min                  | Q23 = 8 L/min
            |                                v
            |                           +----------+
            +-------------------------- |  TANK 3  | <--- Inlet 3 (5 L/min, 20 g/L)
                                        +----------+
                                             |
                                             v Exit 3 (10 L/min)
        """

    def solve(self) -> List[float]:
        return MatrixOps.solve_gaussian(self.M, self.b)


# -----------------------------------------------------------------------------
# INTERACTIVE CLI & AUTO-GRADER SYSTEM
# -----------------------------------------------------------------------------
class EngineeringSimulationRunner:
    """Executes educational task flows, unit tests, and auto-grading routines."""

    @staticmethod
    def run_tests() -> bool:
        """Executes suite of tests for CI/CD integration on GitHub Actions."""
        print("=================================================================")
        print("RUNNING AUTOMATED UNIT TESTS FOR ENGINEERING MATRIX SIMULATION")
        print("=================================================================")
        
        all_passed = True
        
        # Test 1: Mechanical Truss
        truss = MechanicalTrussTask()
        sol_truss = truss.solve()
        expected_f1 = 42.4264  # f1 force component magnitude
        if abs(sol_truss[0] - expected_f1) < 0.1:
            print("[PASS] Task 1: Mechanical Truss Solver (f1 = {:.2f} kN)".format(sol_truss[0]))
        else:
            print("[FAIL] Task 1: Mechanical Truss Solver. Got {:.2f}, Expected {:.2f}".format(sol_truss[0], expected_f1))
            all_passed = False

        # Test 2: Electrical Circuit
        elec = ElectricalCircuitTask()
        sol_elec = elec.solve()
        expected_v1 = 19.0303
        if abs(sol_elec[0] - expected_v1) < 0.1:
            print("[PASS] Task 2: Electrical Circuit Solver (V1 = {:.2f} V)".format(sol_elec[0]))
        else:
            print("[FAIL] Task 2: Electrical Circuit Solver. Got {:.2f}, Expected {:.2f}".format(sol_elec[0], expected_v1))
            all_passed = False

        # Test 3: Fluid Mixing Tanks
        tanks = FluidMixingTanksTask()
        sol_tanks = tanks.solve()
        expected_c1 = 64.1509
        if abs(sol_tanks[0] - expected_c1) < 0.1:
            print("[PASS] Task 3: Fluid Mixing Tanks Solver (c1 = {:.2f} g/L)".format(sol_tanks[0]))
        else:
            print("[FAIL] Task 3: Fluid Mixing Tanks Solver. Got {:.2f}, Expected {:.2f}".format(sol_tanks[0], expected_c1))
            all_passed = False

        print("-----------------------------------------------------------------")
        if all_passed:
            print("ALL AUTOMATED TESTS PASSED SUCCESSFULLY! (CI/CD Exit Code 0)")
        else:
            print("SOME TESTS FAILED! (CI/CD Exit Code 1)")
        print("=================================================================")
        return all_passed

    @staticmethod
    def print_solutions():
        """Prints complete mathematical breakdown for all 3 engineering tasks."""
        print("\n=================================================================")
        print("COMPLETE ENGINEERING MATRIX SOLUTIONS & MATHEMATICAL DERIVATIONS")
        print("=================================================================\n")
        
        # Task 1
        t1 = MechanicalTrussTask()
        print(t1.title)
        print(t1.display_ascii())
        print("Equilibrium System A * x = b:")
        for r in t1.A:
            print("  ", [round(val, 3) for val in r])
        print("Vector b:", t1.b)
        det_A = MatrixOps.determinant(t1.A)
        print(f"Determinant det(A) = {det_A:.4f} (Statically Determinate & Stable)")
        s1 = t1.solve()
        print("\nCalculated Vector x:")
        for name, val in zip(t1.variable_names, s1):
            state = "TENSION" if val > 0.01 else ("COMPRESSION" if val < -0.01 else "REACTION/NEUTRAL")
            print(f"  • {name:<20}: {val:8.3f} kN")

        print("\n" + "-"*65 + "\n")

        # Task 2
        t2 = ElectricalCircuitTask()
        print(t2.title)
        print(t2.display_ascii())
        print("Nodal Conductance Matrix G (Siemens):")
        for r in t2.G:
            print("  ", r)
        print("Current Injection Vector I:", t2.I)
        s2 = t2.solve()
        print("\nCalculated Node Voltages V = G^(-1) * I:")
        for name, val in zip(t2.variable_names, s2):
            print(f"  • {name:<20}: {val:8.3f} Volts")
        
        Z = t2.get_impedance_matrix()
        print("\nNodal Impedance/Resistance Matrix Z = G^(-1) (Ohms):")
        for r in Z:
            print("  ", [round(val, 4) for val in r])

        print("\n" + "-"*65 + "\n")

        # Task 3
        t3 = FluidMixingTanksTask()
        print(t3.title)
        print(t3.display_ascii())
        print("Flow Rate Coefficient Matrix M:")
        for r in t3.M:
            print("  ", r)
        print("Mass Inlet Vector b:", t3.b)
        s3 = t3.solve()
        print("\nCalculated Steady-State Solute Concentrations c:")
        for name, val in zip(t3.variable_names, s3):
            print(f"  • {name:<20}: {val:8.3f} g/L")
        print("\n=================================================================\n")

    @staticmethod
    def run_interactive():
        """Interactive student task mode."""
        print("\n" + "="*65)
        print("  WELCOME TO THE ENGINEERING MATRIX SOLUTIONS INTERACTIVE LAB  ")
        print("="*65)
        print("Select an Engineering Task to solve:\n")
        print("  [1] Task 1: The Mechanical Truss Force Equilibrium")
        print("  [2] Task 2: The Electrical Circuit Nodal Voltage Matrix")
        print("  [3] Task 3: The Fluid Mixing Tanks Mass Balance Matrix")
        print("  [4] Display Full Worked Solutions & Matrix Properties")
        print("  [5] Exit")
        
        try:
            choice = input("\nEnter choice [1-5]: ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\nExiting.")
            return

        if choice == '1':
            t1 = MechanicalTrussTask()
            print("\n" + "="*65)
            print(t1.title)
            print("="*65)
            print(t1.display_ascii())
            print("PROBLEM STATEMENT:")
            print("Calculate the member forces f1, f2, f3 and reaction forces R1x, R1y")
            print("for a applied point load P = 60 kN at Joint 2 with 45-degree angles.")
            print("\nHint: Matrix equation is A * x = b.")
            sol = t1.solve()
            try:
                ans_f1 = float(input("\nEnter your calculated force for Member 1 (f1 in kN): "))
                if abs(ans_f1 - sol[0]) < 1.0:
                    print("CORRECT! f1 = {:.2f} kN".format(sol[0]))
                else:
                    print("INCORRECT. The correct value for f1 is {:.2f} kN.".format(sol[0]))
            except ValueError:
                print("Invalid numerical input.")
            print("\nFull Task 1 Solution Vector:")
            for name, val in zip(t1.variable_names, sol):
                print(f"  {name}: {val:.3f}")

        elif choice == '2':
            t2 = ElectricalCircuitTask()
            print("\n" + "="*65)
            print(t2.title)
            print("="*65)
            print(t2.display_ascii())
            print("PROBLEM STATEMENT:")
            print("Determine node voltages V1, V2, V3 given node current injections")
            print("I = [12A, 0A, 4A]^T and conductance matrix G.")
            sol = t2.solve()
            try:
                ans_v1 = float(input("\nEnter your calculated Node Voltage V1 (in Volts): "))
                if abs(ans_v1 - sol[0]) < 0.5:
                    print("EXCELLENT! V1 = {:.2f} V".format(sol[0]))
                else:
                    print("NOT QUITE. The correct V1 voltage is {:.2f} V.".format(sol[0]))
            except ValueError:
                print("Invalid numerical input.")
            print("\nFull Task 2 Solution Vector:")
            for name, val in zip(t2.variable_names, sol):
                print(f"  {name}: {val:.3f} V")

        elif choice == '3':
            t3 = FluidMixingTanksTask()
            print("\n" + "="*65)
            print(t3.title)
            print("="*65)
            print(t3.display_ascii())
            print("PROBLEM STATEMENT:")
            print("Find the steady-state concentrations c1, c2, c3 in g/L for the 3 CSTR tanks.")
            sol = t3.solve()
            try:
                ans_c1 = float(input("\nEnter your calculated concentration for Tank 1 (c1 in g/L): "))
                if abs(ans_c1 - sol[0]) < 1.0:
                    print("SPOT ON! c1 = {:.2f} g/L".format(sol[0]))
                else:
                    print("INCORRECT. The correct concentration c1 is {:.2f} g/L.".format(sol[0]))
            except ValueError:
                print("Invalid numerical input.")
            print("\nFull Task 3 Solution Vector:")
            for name, val in zip(t3.variable_names, sol):
                print(f"  {name}: {val:.3f} g/L")

        elif choice == '4':
            EngineeringSimulationRunner.print_solutions()

        elif choice == '5':
            print("Goodbye!")


# -----------------------------------------------------------------------------
# MAIN CLI ENTRY POINT
# -----------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Engineering Matrix Solutions Simulation & Student Task Runner"
    )
    parser.add_argument("--test", action="store_true", help="Run automated test suite (for GitHub Actions CI/CD)")
    parser.add_argument("--demo", action="store_true", help="Display full worked engineering solutions")
    parser.add_argument("--interactive", action="store_true", help="Launch interactive CLI task solver")
    parser.add_argument("--task", type=int, choices=[1, 2, 3], help="Execute a specific task directly (1, 2, or 3)")

    args = parser.parse_args()

    if args.test:
        success = EngineeringSimulationRunner.run_tests()
        sys.exit(0 if success else 1)
    elif args.demo:
        EngineeringSimulationRunner.print_solutions()
    elif args.task:
        if args.task == 1:
            t = MechanicalTrussTask()
            print(f"Executing {t.title} Solution: {t.solve()}")
        elif args.task == 2:
            t = ElectricalCircuitTask()
            print(f"Executing {t.title} Solution: {t.solve()}")
        elif args.task == 3:
            t = FluidMixingTanksTask()
            print(f"Executing {t.title} Solution: {t.solve()}")
    elif args.interactive:
        EngineeringSimulationRunner.run_interactive()
    else:
        # Default behavior: Print solutions & hint how to run interactive/test
        print("=================================================================")
        print("   ENGINEERING MATRIX SOLUTIONS SIMULATION PACKAGE (GITHUB READY)")
        print("=================================================================")
        print("Usage Options:")
        print("  python3 matrix_engineering_simulation.py --interactive   Launch Student Interactive Mode")
        print("  python3 matrix_engineering_simulation.py --test          Run GitHub Actions CI/CD Tests")
        print("  python3 matrix_engineering_simulation.py --demo          View Detailed Matrix Solutions")
        print("  python3 matrix_engineering_simulation.py --task [1|2|3]  Run Specific Task Scenario")
        print("=================================================================\n")
        EngineeringSimulationRunner.print_solutions()

if __name__ == "__main__":
    main()
