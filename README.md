# Numerical-PDE-project

###Overview

This project aims to find numerical solutions for a system of 3D partial differential equations defined on the domain Ω = (-1, 1) × (-1, 1) × (-1, 1). The system of equations involves six different fields (F1 to F6), with boundary conditions and parameters set for different values of σ (0.1, 1, 10, 100). The project uses the Symmetric Gauss-Seidel (SGS) method with Newton iteration to solve the system and analyze the convergence of the solution.

### Symmetric Gauss-Seidel (SGS) with Newton Iteration
- Discretization: 100³ mesh points with step size 0.02
- Convergence threshold: 10⁻⁵
- Newton iteration for residual reduction
- Forward and backward scans for comprehensive convergence

### Key Features
- Adaptive iteration process
- Robust convergence for varying σ values
- Efficient handling of 3D spatial discretization
- Dynamic residual monitoring

## Results
### Convergence Analysis
| σ Value | Iterations to Converge |
|---------|----------------------|
| 0.1     | 2                    |
| 1       | 6                    |
| 10      | 88                   |
| 100     | 1319                 |

### Key Findings
- Energy diffusion behavior varies significantly with σ
- Higher σ values show restricted energy spread
- Convergence rate inversely proportional to σ
- Method remains stable across all test cases

## Implementation
### Dependencies
- MATLAB (tested on R2022a)

### Usage
1. Clone the repository
2. Run the main script:
```matlab
[iter_arr, res_arr, A1, A2, A3, A4, A5, A6] = a1_code(Nx, Ny, Nz, sigma, f1, f2, f3, f4, f5, f6, thres)
