# Chapter 6: Numerical Differentiation and ODE

## Numerical Differentiation

### Differentiate interpolating polynomial

Interpolate the function then differentiate the function

!!! note "Error analysis"
    $$
    \begin{align*}
    f'(x_i) - p'(x_i) = \frac{1}{(n+1)!}f^{(n+1)}(\xi_x)\prod^n_{j=0, j \neq i}(x_i-x_j)\\
    \xi \in (a, b)
    \end{align*}
    $$

    $$
    \begin{align*}
    f'(x) - p'(x) = \frac{1}{n!}f^{(n+1)}(\eta)\prod^{n-1}_{i=0}(x-\xi_i)\\
    \xi_i \in (x_i, x_{i+1}), \eta \text{ in interval containing } \{ \xi_0, \cdots, \xi_{n-1}, x \}
    \end{align*}
    $$

    ??? "Raw image ptsd"

        ![ch6 polynomial error 1](../images/ch6_polynomial_error_1.png)
     
        ![ch6 polynomial error 2](../images/ch6_polynomial_error_2.png)

### Finite difference schemes

$\boxed{f'(x) \approx \frac{f(x+h)-f(x)}{h}}$ is order 1 accuracy, because $f'(x) = \frac{f(x+h)-f(x)}{h} - \frac{h}{2}f''(\xi)$

$\frac{h}{2}f''(\xi)$ can give us a bound on the error.

!!! warning "Subtractive cancellation"
    If we keep on decreasing $h$, we will not neccessarily get better accuracy due to **subtractive cancellation**.
    
    The number of significant digits drops as difference of very close numbers are being calculated.
    
    **Roundoff error** prevents us from getting better approximations by decreasing $h$.

$\boxed{f'(x) \approx \frac{f(x+h)-f(x-h)}{2h}}$ is order 2 accuracy

??? question "Proof"
    ![ch6_central_difference_proof.png](../images/ch6_central_difference_proof.png)

### Richardson Extrapolation

Let $\varphi(h) = \frac{1}{2h}(f(x+h)-f(x-h))$

Richardson Extrapolation 4th order approximation: $\boxed{\frac{4}{3}\varphi(h/2)-\frac{1}{3}\varphi(h)}$

??? question "Proof"
    ![ch6_richardson_proof.png](../images/ch6_richardson_proof.png)

## IVP (Initial Value Problem)

Given that $\boxed{\varphi'(t) = x' = f(t, x), x(t_0) = x_0}$, obtain a solution $\boxed{x = \varphi(t)}$

By numerical solutions, we mean a sequence of approximated values of $x$ at a sequence of values of $t$.

Usually, denote $h = t_i - t_{i-1}$

### Euler's method

$$
\varphi(t) \approx x_t = x_0 + f(t_0, x_0)(t - t_0)
$$

??? example
    ![ch6_euler_method_eg.png](../images/ch6_euler_method_eg.png)

Global truncation error is $O(h)$. Because the local truncation error in each step is $\frac{1}{2}f''(\xi_i)h^2$, and we take $(t_i-t_0)/h$ steps.

Can achieve higher order, but it's complicated to derive formula for higher derivatives.

### Taylor Series method

Expand the taylor series of $\varphi(t = t_0 + h)$, then express the derivatives in terms of lower derivatives.

!!! example
    ![ch6_taylor_series.png](../images/ch6_taylor_series.png)

### Runge-Kutta method (2nd order)

$$
\begin{align*}
K_1 &= f(t_i, x_i)\\
\bar{x}_{i+1} &= x_i + hf(t_i, x_i)\\
K_2 &= f(t_i + h, x_i + hK_1)\\
x_{i+1} &= x_i + h(K_1/2 + K_2/2)
\end{align*}
$$

??? example
    ![ch6_runge_kutta_eg.png](../images/ch6_runge_kutta_eg.png)

??? "Derivation (trapezoidal rule)"
    ![ch6_runge_kutta_derive.png](../images/ch6_runge_kutta_derive.png)

Local truncation error: $O(h^3)$

??? Question "Proof"
    ![ch6_runge_kutta_error.png](../images/ch6_runge_kutta_error.png)

In general, can choose

$$
\begin{align*}
K_1 &= f(t_i, x_i)\\
\bar{x}_{i+1} &= x_i + bhf(t_i, x_i)\\
K_2 &= f(t_i + ah, x_i + bhK_1)\\
x_{i+1} &= x_i + h(c_1K_1 + c_2K_2)\\
c_1 + c_2 = 1, & ac_2 = 1/2 = bc_2
\end{align*}
$$

??? "Derivation (Taylor series)"
    ![ch6_runge_kutta_general_derive.png](../images/ch6_runge_kutta_general_derive.png)

## BVP (Boundary Value Problem)

Given that $\boxed{\varphi''(t) = x'' = f(t, x, x'), x(0) = \alpha, x_{n+1} = \beta}$, obtain a solution $\boxed{x = \varphi(t)}$

Approximate

- $\boxed{f'(x) \approx \frac{1}{2h} (f(x+h)-f(x-h))}$
- $\boxed{f''(x) \approx \frac{1}{h^2} (f(x+h)-2f(x)+f(x-h))}$, order 2

!!! important "Discretized BVP"
    ![ch6_bvp_1.png](../images/ch6_bvp_1.png)

Linear case: $\boxed{x'' = f(t, x, x') = u(t) + v(t)x + w(t)x'}$

!!! important "Discretized BVP (Linear case"
    ![ch6_bvp_2.png](../images/ch6_bvp_2.png)

Coefficient matrix is a tridiagonal matrix. It is strictly diagonally dominant (proof omitted) so inverse exists. It can be solved in $O(n)$.

## Convergence analysis

Basically, both Euler method and BVP converge when $h \rightarrow 0$ and $n \rightarrow 0$. But we don't have to prove it.

??? abstract "Convergene of Euler method"
    ![ch6_convergence_euler.png](../images/ch6_convergence_euler.png)

??? abstract "Convergene of finite difference BVP"
    ![ch6_convergence_bvp.png](../images/ch6_convergence_bvp.png)
