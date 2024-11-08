# Numerical-PDE-project

**Overview**

This project aims to find numerical solutions for a system of 3D partial differential equations defined on the domain Ω = (-1, 1) × (-1, 1) × (-1, 1). The system of equations involves six different fields (F1 to F6), with boundary conditions and parameters set for different values of σ (0.1, 1, 10, 100). The project uses the Symmetric Gauss-Seidel (SGS) method with Newton iteration to solve the system and analyze the convergence of the solution.

**Numerical Methods**

To solve the system, the following numerical method is used:
Symmetric Gauss-Seidel (SGS) Method with Newton Iteration: This method provides fast convergence, making it suitable for large values of σ. The SGS method is applied iteratively, with Newton's method used to reduce residuals and ensure convergence for all values of σ. The convergence criteria are set based on a tolerance value of 10^-5.

**Results**

The numerical results are presented for different values of σ:
- Figures: The report includes figures showing the distribution of field values for each value of σ. As σ increases, it becomes harder for the non-zero values to spread across the domain.
- Convergence Analysis: A summary table and plots showing the number of iterations needed for convergence for each value of σ. The method converges faster for smaller values of σ, while larger values require more iterations.

