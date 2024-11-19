# Numerical-PDE-project

### Overview

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
```

### Example code
```matlab
% Set parameters
Nx = 100; Ny = 100; Nz = 100;  % Grid points
sigma = 1;                      % Diffusion parameter
thres = 1e-5;                  % Convergence threshold

% Define boundary conditions
f1 = @(y,z) double(abs(y)<=0.2 & abs(z)<=0.2);  % F₁(-1,y,z)
f2 = @(y,z) zeros(size(y));                      % F₂(1,y,z)
f3 = @(x,z) double(abs(x)<=0.2 & abs(z)<=0.2);  % F₃(x,-1,z)
f4 = @(x,z) zeros(size(x));                      % F₄(x,1,z)
f5 = @(x,y) double(abs(x)<=0.2 & abs(y)<=0.2);  % F₅(x,y,-1)
f6 = @(x,y) zeros(size(x));                      % F₆(x,y,1)

% Run solver
[iter_arr, res_arr, A1, A2, A3, A4, A5, A6] = a1_code(Nx, Ny, Nz, sigma, ...
    f1, f2, f3, f4, f5, f6, thres);

% Display results
fprintf('Converged in %d iterations\n', length(iter_arr));
fprintf('Final residual: %e\n', res_arr(end));
