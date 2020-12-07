# Chapter 4: Approximating functions

## Function interpolation

!!! question "Motivating question"
    Given a set of points $(x_i, y_i)$ for $i = 0, 1, \cdots, n$, where the $x_i$ are **distinct** values of the independent variable and the $y_i$ are the corresponding value of some function $f$ at $x_i$, that is $y_i = f(x_i)$.

    - Determine a function $g$ that satisfies $g(x_i) = f(x_i)$ and predicts the value of $f$ at some value $x$ not listed among $x_i$

!!! success "Goal"
    Interpolation: The function $g$ is determined by requiring the **error to be zero at each of the $x_i$**, i.e. $g(x_i) = y_i$.

Linear interpolation: $g(x; a_0, \cdots, a_n) = a_0g_0(x) + a_1g_1(x) + \cdots + a_ng_n(x)$

Typical linear interpolation methods include:

- Polynomial interpolation
- Spline (piecewise polynomial)
- ~~Trignometric interpolation~~

### Polynomial interpolation

Goal: Seek a polynomial $p$ of the lowest degree for which $p(x_i) = y_i, (0 \leq i \leq n)$

- Polynomials are smooth functions
- Easy to store, manipulate, take derivative/antiderivative

However, polynomials are **rigid** and may create problems at boundaries (Re: Runge's phenomenon)

!!! important "Theorem (degree of polynomial interpolation)"
    If $x_0, x_1, \cdots, x_n$ are **distinct** real numbers, then for arbitrary values $y_0, y_1, \cdots, y_n$, there is a **unique** polynomial $p_n$ of degree **at most $n$** such that $p_n(x_i) = y_i$.

??? question "Proof (uniqueness)"
    Suppose there were two such polynomails $p_n$ and $q_n$. Then $r_n(x) = (p_n - q_n)(x_i) = 0, 0 \leq i \leq n$.

    Since the degree of $r_n$ can be at most $n$, this polynomial can have at most $n$ zeros if it is not the $0$ polynomial. Since the $x_i$ are distinct, $r_n$ has $n+1$ zeros; it must therefore be $0$. Hence $p_n \equiv q_n$.

!!! question "Proof (existence)"
    The proof for existence actually let us know some ways to generate the polynomial. So we leave it to the next part.

!!! important "Theorem (Polynomial interpolation error)"
    To each $x$ there corresponds a point $\eta_x \in (a, b)$ such that:
    
    $$
    f(x) - p(x) = \frac{1}{(n+1)!} f^{(n+1)}(\eta_x) \prod^n_{i=0}(x-x_i)
    $$

??? question "Proof"
    $x$ is fixed.

    Let $w(t) = \prod^n_{i=0} (t-x_i)$, $\phi(t) \equiv f(t) - p(t) + \lambda w(t)$, where $\phi(x) = 0$ and $\lambda = \frac{f(x)-p(x)}{w(x)}$

    $\phi \in C^{n+1}[a, b]$, and $\phi$ vanishes at the $n+2$ points $x, x_0, x_1, \cdots, x_n$. 

    By Mean-Value Theorem, $\phi'$ has at least $n+1$ distinct zeros in $[a, b]$ (one zero in $\phi'$ between each original zero in $phi$)

    Similarly, $\phi''$ has at least $n$ distinct zeros in (a, b).

    Eventually, $\phi^{(n+1)}$ has at least one zero, say $\eta_x$ in $(a, b)$.

    Now $w^{(n+1)}(t)$ = \frac{d^{n+1}}{dt^{n+1}}\prod^n_{i=1}(t-x_i) = (n+1)!$ (polynomial in $t$)

    $\phi^{(n+1)} = f^{(n+1)} - p^{(n+1)} - \lambda w^{(n+1)} = f^{(n+1)} - (n+1)!\lambda$ ($p$ has degree at most $n$)

    Lastly, $0 = \phi^{(n+1)}(\eta_x) = f^{(n+1)}(\eta_x) - (n+1)!\frac{f(x)-p(x)}{w(x)}$


!!! note "Wolfram alpha command"
    ```
    interpolate [(5,1),(-7,-23),(-6,-54), (0, -954)]
    ```

### Lagrange interpolation

$$
\begin{align*}
p_n(x) = y_0 \ell_0(x) + y_1 \ell_1(x) + \cdots + y_n \ell_n(x) \\
\ell_i(x) = \prod^n_{j=0, j \neq i} \frac{x-x_j}{x_i-x_j}
\end{align*}
$$

where $\ell_i$ are polynomials that depend on $x_0, \cdots x_n$ but not on $y_0 \cdots y_n$.

Note that $y_i = p_n(x_i) = y_0 \ell_0(x_i) + y_1 \ell_1(x_i) + \cdots + y_n \ell_n(x_i)$.

So $\ell_i(x_j) = 1$ if $i = j$ and $\ell_i(x_j) = 0$ if $i \neq j$

!!! example
    $(1, 1), (2, 8), (3, 27)$

    The interpolating polynomial is:

    $1 \cdot \frac{x-2}{1-2} \cdot \frac{x-3}{1-3} + 8 \cdot \frac{x-1}{2-1} \cdot \frac{x-3}{2-3} + 27 \cdot \frac{x-1}{3-1} \cdot \frac{x-2}{3-2} = 6x^2 - 11x + 6$

    ```
    interpolate [(1, 1) (2, 8) (3, 27)]
    ```

!!! note
    $\sum^n_{i=0}\ell_i(x) = 1$

    Try interpolating $P(x) = f(x) = 1$

While elegant and theoretically important, the Lagrange interpolation formula is not suitable for actual computations, because:

- the evaluation of the Lagrange interpolation requires more floating point operations than necessary (Re: Horner's method)
- difficult to obtain derivative
- the inlcusion of additional data points require the **reevalutaion** of all Lagrange polynomials
- **Runge's phenomenon**: If the points are equally spaced, then there is large oscillation above and below the true function.

### Newton interpolation

$$
\begin{align*}
p_0(x_0) = c_0 = y_0 \\
p_k(x) = p_{k-1}(x) + c_k(x-x_0)(x-x_1)\cdots(x-x_{k-1})\\
c_k = \frac{y_k - p_{k-1}(x)}{(x_k-x_0)(x_k-x_1)\cdots(x_k-x_{k-1})}
\end{align*}
$$

Divided differences notation

$$
\begin{align*}
c_0 = p_n(x_0) = y_0 = [y_0]\\
c_1 = \frac{y_1 - y_0}{x_1 - x_0} = \frac{[y_1] - [y_0]}{x_1 - x_0} = [y_0, y_1]\\
c_k = \frac{[y_1, \cdots, y_k] - [y_0, \cdots, y_{k-1}]}{x_k - x_0} = [y_0, \cdots, y_k]\\
p_n(x) = [y_0] + \sum^n_{k=1}[y_0, y_1, \cdots, y_k](x-x_0) \cdots (x-x_{k-1})
\end{align*}
$$

!!! example

    $(0, 1), (1, 3), (3, 2)$

    |            | $k=0$  | $k=1$  | $k=2$  |
    |------------|--------|--------|--------|
    | $x_0 = 0$  | $[y_0] = 1$*  |   |   |
    | $x_1 = 1$  | $[y_1] = 3$  | $[y_0, y_1] = 2$*  |   |
    | $x_2 = 3$  | $[y_2] = 2$  | $[y_1, y_2] = -1/2$  | $[y_0, y_1, y_2] = -5/6$*  |

    $p_n = [y_0] + [y_0, y_1](x-x_0) + [y_0, y_1, y_2](x-x_0)(x-x_1)$

    $p_n = 1 + 2(x-x_0) + -\frac{5}{6}(x-x_0)(x-x_1) = -\frac{5}{6}x^2 + \frac{17}{6}x + 1$

### Hermite interpolation

(skipped)

### Spline interpolation

A spline function consists of polynomial pieces on subintervals joined together with certain continuity conditions.

Formally suppose that $n+1$ points $t_0, t_1, \cdots, c_n$ have been specified and satisfy $t_0 < t1 < \cdots < t_n$. These points are called **knots**.

Suppose also that an integer $k \geq 0$ has been perscribed. A spline function of degree $k$ having knots $t_0, \cdots, t_n$ is a function $S$ such that:

- On each interval $[t_{i-1}, t_i]$, $S$ is a polynomial of degree $\leq k$
- $S$ has a continuous $(k-1)$st derivative on $[t_0, t_n]$

i.e. Piecewise polynomial of degree at most $k$ having continuous derivatives of all orders up to $k-1$.

- Spline of degree 0: Piecewise constant, step funtion (right continuous), intervals do not intersect
- Spline of degree 1: Piecewise linear, continuous

Checking spline: Left hand and right hand limit of value and (deg-1) derivatives are continuous

#### Cubic spline

!!! abstract "Derivation"
    $z_i = S''(t_i)$

    $S_i^{''}(x) = A_ix+B_i$

    $A_it_i + B_i = z_i, A_it_{i+1} + B_i = z_{i+1}$

    $h_i = t_{i+1} - t_i$

    $A_i = \frac{z_{i+1}-z_i}{h_i}, B_i = z_i - \frac{z_{i+1} - z_i}{h_i} t_i$

    $S_i^{''}(x) = A_ix+B_i = \frac{z_{i+1}}{h_i}(x-t_i) + \frac{z_i}{h_i}(t_{i+1}-x)$, straight line between $z_i$ and $z_{i+1}$

    $S_i^{'}(x) = -\frac{1}{2}\frac{z_i}{h_i}(t_{i+1}-x)^2 + \frac{1}{2}\frac{z_{i+1}}{h_i}(x-t_i)^2+C_i$ (integration)

    $S_i(x) = \frac{z_i}{6h_i}(t_{i+1}-x)^3 + \frac{z_{i+1}}{6h_i}(x-t_i)^3 + C_i(x-t_i) + D_i$

    Interpolation conditions $S_i(t_i) = y_i$ and $S_i(t_{i+1}) = y_{i+1}$, find $C_i$ and $D_i$: $D_i = y_i - \frac{z_ih_i^2}{6}, C_i = \frac{y_{i+1}-y_i}{h_i} + \frac{h_i}{6}(z_i-z_{i+1})$

    **Determining $z_1, z_2, \cdots, z_{n-1}$:**

    Use continuity conditions for $S'$ at the interior knots $t_i$: $S_{i-1}^{'}(t_i) = S_{i}^{'}(t_i)$

    TODO

!!! important "Theorem (Optimality of Natural Cubic Spline)"
    Let $f^{''}$ be continuous on $[a, b]$ and let $a = t_0 < t_1 < \cdots < t_n = b$. If $S$ is the natural cubic spline interpolating $f$ at the knots $t_i$, then $\int^b_a{[S^{''}(x)]^2} dx \leq \int^b_a{[f^{''}(x)]^2} dx$

??? question "Proof"
    Let $g \equiv f-S$.

    $(f^{''})^2 = (S^{''} + g^{''})^2 = (S^{''})^2 + (g^{''})^2 + 2S^{''}g^{''}$

    $\int^a_b{(f^{''})^2}dx = \int^a_b{(S^{''})^2}dx + \int^a_b{(g^{''})^2}dx + \int^a_b{2S^{''}g^{''}}dx$

    Note that $\int^a_b{(g^{''})^2}dx \geq 0$. So to prove the inequality, we can choose to show that $\int^a_b{2S^{''}g^{''}}dx \geq 0$.
    
    Idea: Integration by parts

    $\int^a_b{2S^{''}g^{''}}dx = \sum^n_{i=1}{\int^{t_{i}}_{t_{i-1}}{2S^{''}g^{''}}dx} = \sum^n_{i=1}{S^{''}g^{'}(t_i) - {S^{''}g^{'}(t_{i-1}) - \int^{t_{i}}_{t_{i-1}}{2S^{'''}g^{'}}dx}$

    TODO

## Function approximation

### Least squares

$$
\begin{align*}
X\beta = Y \\
\hat{\beta} = (X^TX)^{-1}X^Ty
\end{align*}
$$

### Function approximation over interval

!!! example
    Let $f(x) = e^x$, let $p(x) = \alpha_0 + \alpha_1x$. Approximate $f(x)$ over $[-1, 1]$. 

    So choose $\alpha_0, \alpha_1$ to minimize $g(\alpha_0, \alpha_1) = \int_{-1}^1{[e^x - \alpha_0 - \alpha_1]^2 dx}$

    $0 = \frac{\partial g}{\partial \alpha_0} = 2 \int_{-1}^1{[e^x - \alpha_0 - \alpha_1] (-1)} dx$

    $0 = \frac{\partial g}{\partial \alpha_0} = 2 \int_{-x}^1{[e^x - \alpha_0 - \alpha_1] (-1)} dx$

    Solving gives $2\alpha_0 = e-e^{-1}$, $\frac{2}{3}\alpha_1 = 2e^{-1}$.

#### Gram-Schmidt process

Something related to orthogonal functions...

E.g. Orthogona function: Fourier series

TODO