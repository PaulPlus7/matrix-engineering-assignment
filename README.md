# Engineering Matrix Solutions — Interactive Simulation & Lab Assignment

[![Autograding Tests](../../workflows/Autograding%20Tests/badge.svg)](../../actions)
**Level 4 HNC / Undergraduate Engineering Mathematics**  
*Linear Algebra & Matrix Applications in Engineering*

---

## 📌 Assignment Overview

Welcome to the **Engineering Matrix Solutions** lab assignment! In this practical module, you will apply fundamental linear algebra concepts—specifically matrix inversion, matrix multiplication, and system equilibrium solving ($A \mathbf{x} = \mathbf{b}$)—to solve real-world problems in structural, electrical, and chemical engineering.

This repository contains an interactive Python simulation environment (`matrix_engineering_simulation.py`) that models three real engineering systems:

1. **The Mechanical Truss** (Structural Engineering — Force Equilibrium)
2. **The Electrical Circuit** (Electrical Engineering — Nodal Conductance Matrix & KCL)
3. **The Fluid Mixing Tanks** (Chemical/Process Engineering — Continuous Mass Balance)

---

## 🛠️ Getting Started & Requirements

### Prerequisites
* **Python 3.8 or higher** (Standard Library only — no external dependencies required!)
* **Git** installed on your machine
* An active **GitHub** account linked to GitHub Classroom

### Running the Simulation

Clone this repository to your local computer or open it in **GitHub Codespaces**:

```bash
git clone https://github.com/YourOrganization/engineering-matrix-solutions.x.git
cd engineering-matrix-solutions
```

You can execute the script using several flags depending on what you want to do:

| Command Flag | Description |
| :--- | :--- |
| `python3 matrix_engineering_simulation.py --interactive` | **Student Interactive Mode** — Guided menu with ASCII schematics and problem-solving prompts. |
| `python3 matrix_engineering_simulation.py --test` | **Automated Test Runner** — Runs unit tests locally to verify your answers before submitting. |
| `python3 matrix_engineering_simulation.py --demo` | **Full Worked Solutions** — Displays step-by-step mathematical derivations and full solution matrices. |
| `python3 matrix_engineering_simulation.py --task 1` | Runs Task 1 (Mechanical Truss) individually. |
| `python3 matrix_engineering_simulation.py --task 2` | Runs Task 2 (Electrical Circuit) individually. |
| `python3 matrix_engineering_simulation.py --task 3` | Runs Task 3 (Fluid Mixing Tanks) individually. |

---

## 📐 Lab Tasks Breakdown

### 🏗️ Task 1: The Mechanical Truss (Structural Engineering)
* **Problem:** Analyze a 2D pin-jointed bridge truss subjected to a downward point load $P = 60\text{ kN}$ at Joint 2.
* **Governing Equation:** $A \mathbf{x} = \mathbf{b}$, where $A$ is a $5 \times 5$ equilibrium matrix of direction cosines, $\mathbf{x} = [f_1, f_2, f_3, R_{1x}, R_{1y}]^T$, and $\mathbf{b} = [0, 0, 0, -60, 0]^T$.
* **Key Questions:**
  1. Calculate internal member forces $f_1, f_2, f_3$ and reaction forces $R_{1x}, R_{1y}$.
  2. Verify if the truss is statically determinate by evaluating $\det(A)$.
  3. Classify internal member forces as **Tension** ($f > 0$) or **Compression** ($f < 0$).

### ⚡ Task 2: The Electrical Circuit (Electrical Engineering)
* **Problem:** Determine nodal voltages in a 3-node resistor network powered by DC current injections ($I_1 = 12\text{ A}, I_3 = 4\text{ A}$).
* **Governing Equation:** $G \mathbf{v} = \mathbf{i}$, where $G$ is the $3 \times 3$ nodal conductance matrix ($G_{ij} = -1/R_{ij}$, $G_{ii} = \sum 1/R$).
* **Key Questions:**
  1. Construct matrix $G$ using branch conductances (Siemens).
  2. Solve for node voltages $\mathbf{v} = G^{-1} \mathbf{i}$.
  3. Compute the nodal impedance matrix $Z = G^{-1}$ and state which node exhibits the highest equivalent resistance to ground.

### 🧪 Task 3: The Fluid Mixing Tanks (Chemical Engineering)
* **Problem:** Model steady-state solute concentrations in a Continuous Stirred-Tank Reactor (CSTR) system consisting of 3 interconnected mixing tanks with recirculation loops.
* **Governing Equation:** $M \mathbf{c} = \mathbf{b}$, representing continuous mass balance ($\text{Mass In} = \text{Mass Out}$).
* **Key Questions:**
  1. Formulate the flow rate coefficient matrix $M$.
  2. Calculate steady-state concentrations $c_1, c_2, c_3$ in g/L.
  3. Explain why Tanks 1 and 2 achieve identical steady-state concentrations.

---

## 🤖 Automated Grading & Continuous Integration (CI/CD)

This repository is integrated with **GitHub Actions**. Every time you push code or commits to GitHub:

1. GitHub automatically executes the workflow defined in `.github/workflows/grader.yml`.
2. The autograder runs `python3 matrix_engineering_simulation.py --test`.
3. You will receive immediate feedback directly on GitHub (a green tick `✔` for pass or a red cross `✖` if solutions deviate from expected engineering tolerances).

---

## 📤 Submission Instructions

1. Complete all calculations using interactive mode or your own Python matrix solver script.
2. Run `python3 matrix_engineering_simulation.py --test` locally to verify that all test suites pass.
3. Commit and push your work to GitHub:
   ```bash
   git add .
   git commit -m "Submit completed engineering matrix solutions"
   git push origin main
   ```
4. Check the **Actions** tab on your GitHub repository page to verify your test results!
