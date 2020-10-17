# Chapter 3: Numerical Linear Algebra

$$
Ax = b
$$

Easy to solve for lower triangular and upper triangular matrix:

- Lower triangular matrix $Lx = b$: Forward substitution
- Upper triangular matrix $Ux = b$: Backward substitution
- Computational complexity: $O(n^2)$

Reduce general matrix to triangular matrix via **Gaussian Elimination**, complexity $O(n^3)$.

## Gaussian Elimination with Partial Pivoting

- In Gaussian elimination, a row interchange was needed when one of the pivot elements $a^{(k)}_{kk}$ is $0$. 
- When $a^{(k)}_{kk}$ is nonzero but very small, pivoting using this row may cause **round-off error** due to large multipliers.

!!! important "Partial pivoting strategy"
    - Select an element in the same column that is below the diagonal and has the largest absolue value, and use that row as the pivot. Swap rows if neccessary.

## LU Factorization

1. Decompose $A = LU$ ($O(\frac{2}{3}n^3)$, smaller constant)
2. Solve $Ly = b$ ($O(n^2)$)
3. Solve $Ux = y$ ($O(n^2)$)

Finding the decomposition: For each cell in $A$, just solve the multiplication equation. (Re: Doolittle algorithm)

## Iterative methods

### Norm

- $L_2$ norm: $||x||_2 = \sqrt{x_1^2 + x_2^2 + \cdots + x_n^2}$
- $L_1$ norm: $||x||_1 = |x_1| + |x_2| + \cdots + |x_n|$
- $L_\infty$ norm: $||x||_\infty = \max{(|x_1|, |x_2|, \cdots, |x_n|)}$

Matrix norm induced by a given vector norm:

- $||A|| = \sup_{x \neq 0} \frac{||Ax||}{||x||} \implies ||Ax|| \leq ||A||||x||$
- Special cases refer to slides
    - 1-norm is max col sum of abs elements
    - inf-norm is max row sum of abs elements
    - 2-norm is sqrt of largest eigenvalue of $A^TA$

Bound on relative error: $\frac{||\delta x||}{||x||} \leq ||A||||A^{-1}||\frac{||\delta b||}{||b||}$ (smaller is better)

Condition number: $||A||||A^{-1}||$ (smaller is better)

??? Derivation

### Jacobi method

$$
\begin{align*}
A = D + R\\
R = \begin{bmatrix}
0 & a_{12} & a_{13} \\
a_{21} & 0 & a_{23} \\
a_{31} & a_{32} & 0 \\
\end{bmatrix}\\
D = \begin{bmatrix}
a_{11} & 0 & 0 \\
0 & a_{22} & 0 \\
0 & 0 & a_{33} \\
\end{bmatrix}\\
x^{(k+1)} = D^{-1}(b - Rx^{(k)}) = D^{-1} - D^{-1}Rx^{(k)}\\
x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j \neq i}{a_{ij}x_j^{(k)}}\right)
\end{align*}
$$

- Computation of $x_i^{(k+1)}$ requires each element in $x^{(k)}$ except itself.
- The iterative scheme is easy to be parallelized.

!!! important "General convergence condition"
    The standard convergence condition **(for any iterative method)** is when the norm of the iteration matrix (i.e. $||D^{-1}R||$) is less than $1$.

    ??? question "Proof"
        $e^{(k+1)} = X^{(k+1)} - X^* = D^{-1}b - D^{-1}Rx^{(k)}-X^*$

        $(D+R)X^* = b \implies X^* = D^{-1}b - D^{-1}RX^*$

        $e^{(k+1)} = X^{(k+1)} - X^* = D^{-1}R\left(X^* - X^{(k)}\right) = D^{-1}R \cdot e^{(k)} \implies e^{(k+1)} = (D^{-1}R)^{k+1} \cdot e^{(0)}$

!!! important "Strictly diagonally dominant case"
    If the matrix $A$ is **strictly diagonally dominant (i.e. $|a_{ii}| > \sum_{j \neq i}|a_{ij}|)$** (for all rows the absolute value of the diagonal element in a row is strictly greater than than the sum of absolute value of the rest of the elements in that row.), the iterative method is guaraneed to converge.

The Jacobi method sometimes converges even if these conditions are not satisfied.

### Gauss-Seidel method

$$
\begin{align*}
A = L_* + U\\
x^{(k+1)} = L_*^{-1}(b - Ux^{(k)})\\
x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j < i}{a_{ij}x_j^{(k+1)}} - \sum_{j > i}{a_{ij}x_j^{(k)}}\right)
\end{align*}
$$

$L*$ includes diagonal, $U$ does not.

- The computation of $x_i^{(k+1)}$ **uses only the elements of $x^{(k+1)}$ that have already been computed**, and only the elements of $x^{(k)}$ that have not yet advanced to iteration $k+1$.
- Cannot be done in parallel, the values at each iteration are dependent on the order of the original equations.
- Convergence speed faster than Jacobi method.

!!! important "General convergence condition"

    The procedure is known to converge if $A$ is symmetric positive-denite, or $A$ is strictly diagonally dominant.

The Gauss-Seidel method sometimes converges even if these conditions are not satisfied.

### SOR method (not required)

## Gradient descent

Assumptions: Symmetric positive definite

Solving $Ax = b$ is equivalent to minimizing $f(x) = \frac{1}{2}x^TAx-b^Tx+c$

!!! Algorithm

    1. Compute the search direction $d_k = -\nabla f(x_{k-1}) = b - Ax_{k-1}$
    2. Define $x_k = x_{k-1} + \alpha_kd_k$, where $\alpha_k$ is chosen such that $f(x_k) < f(x_{k-1})$
    3. Convergence check. Stop or repeat.

Deciding on $\alpha_k$ is minimizing $f(x_{k-1} + \alpha_kd_k)$ as a function of a single scalar variable $\alpha_k$.

After some derivation of differentiating scalar function $g(\alpha) = f(x_{k-1} + \alpha d_k)$, we obtain $\frac{dg(\alpha)}{d\alpha} = \alpha d_k^TAd_k + d_k^T(A x_{k-1}-b)$

Solving $\frac{dg(\alpha)}{d\alpha} = 0$, obtain $\alpha = \frac{d_k^T(b-Ax_{k-1})}{d_k^TAd_k}$

For gradient descent, choose $d_k = b - Ax_{k-1}$, so $\alpha_k = \frac{||d_k||_2^2}{d_k^TAd_k}$

!!! Algorithm

    1. Compute the search direction $d_k = -\nabla f(x_{k-1}) = b - Ax_{k-1}$
    2. Define $x_k = x_{k-1} + \alpha_kd_k$, where $\alpha_k = \frac{||d_k||_2^2}{d_k^TAd_k}$
    3. Convergence check. Stop or repeat.

Some notes:

- Successive search directions $d_{k+1}$ and $d_{k}$ are orthogonal

---

$R_k = 1-\frac{\sigma_{min}}{\sigma_{max}}$ (smaller is better)

- Smallest and greatest eigenvalue
- Smaller is better (clustered, like circle)
- Large: Eigenvalue scattered, long and thin landscape, zig-zagging on energy landscape
