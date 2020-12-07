# Chapter 2: Solutions of Nonlinear Equations

## Bisection method

- Based on the **intermediate value theorem**. If $f$ is continuous on $[a, b]$, and if $f(a) \cdot f(b) < 0$ (sign flip), then $f$ must have a zero in $(a, b)$. 

### Main loop

1. Set $c = (a+b)/2$ and test whether $f(a)f(c) < 0$.
2. If this is true, then $f$ has a zero in $(a, c)$ (the signs of $f(a)$ and $f(c)$ are different). So we rename $c$ as $b$ and start again with the new interval in $[a, b]$.
3. Else if $f(a)f(c) > 0$, then $f(c)f(b) < 0$ (the sign of $f(a)$ and $f(c)$ are the same). In this case, rename $c$ as $a$.

Repeat until the **length of the interval** or **$f(c)$ is small enough**. 

!!! example
    - Find root of $e^x = sin x$ in the interval $[-4, -3]$, with termination scalar $\epsilon$.
    - Then let $f(x) = e^x - sin x$. Finding zero of $f(x)$ is equivalent to finding root in original equation.
    - $f(-4)f(-3) < 0$
    - Then continue iterating until $f(c) < \epsilon$

!!! important "Therorem (Convergence of Bisection Method)"
    1. The limits $\lim_{n\to\infty} a_n$ and $\lim_{n\to\infty} b_n$ exist.

    ??? question "Proof"
        - $a_0 \leq a_1 \leq \cdots \leq b_0$
        - The sequence $a_n$ is **nondecreasing** and **bounded above**, so $\lim_{n\to\infty} a_n$ exists. Similarly, $\lim_{n\to\infty} b_n$ exists.
    
    2. If $c_n = (a_n + b_n) / 2$ and $\lim_{n\to\infty} c_n = r$, then $|r - c_n| \leq 2^{-(n+1)}(b_0 - a_0)$

    ??? question "Proof"
        - $b_1 - a_1 = \frac{1}{2}(b_0 - a_0)$
        - $b_2 - a_2 = \frac{1}{2}(b_1 - a_1) = \frac{1}{4}(b_0 - a_0)$
        - $\cdots$
        - $b_n - a_n = \frac{1}{2^n}(b_0 - a_0) = 2^{-n}(b_0 - a_0)$

        Thus, $\lim_{n\to\infty}b_n - \lim_{n\to\infty}a_n = \lim_{n\to\infty}2^{-n}(b_0-a_0) = 0$

        Therefore $\lim_{n\to\infty}a_n = \lim_{n\to\infty}b_n \equiv r$

        Since $f$ is continuous, $\lim_{n\to\infty}f(a_n) = f(\lim_{n\to\infty}f(a_n)) = r$ and $\lim_{n\to\infty}f(b_n) = f(\lim_{n\to\infty}f(b_n)) = r$

        $0 \geq f(a_n)f(b_n) \implies 0 \geq [f(r)]^2 \implies f(r) = 0$

        $|r - c_n| \leq \frac{1}{2}(b_n - a_n) = 2^{-(n+1)}(b_0 - a_0)$

### Analysis

- For an error $\epsilon > 0$, if $\frac{1}{2}(b_n - a_n) = 2^{-(n+1)}(b_0 - a_0) \le \epsilon$, we need to choose $n > \frac{\ln[(b_0 - a_0)/\epsilon]}{\ln 2} - 1$

??? Derivation
    - $2^{-(n+1)}(b_0 - a_0) \leq \epsilon$
    - $-(n+1)\ln 2 + \ln(b_0 - a_0) \leq \ln \epsilon$
    - $-\ln 2 + \ln(b_0 - a_0) - \ln \epsilon \leq n \ln 2$
    - $n > \frac{\ln[(b_0 - a_0)/\epsilon]}{\ln 2} - 1$

- Linear convergence: $|r - c_n| \leq \frac{1}{2}(b_n - a_n), |r - c_{n+1}| \leq \frac{1}{2}(\frac{1}{2}(b_n - a_n)) \implies |r - c_{n+1}| \approx \frac{1}{2}|r - c_{n}|$.

- Reliable but **slow**.
- Guaranteed to work if $f(x)$ is continuous.

## Newton-Raphson's Method

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
$$

??? Derivation
    By Taylor's theorem, $0 = f(r) = f(x + h) = f(x) + hf'(x) + \frac{h^2}{2}f''(\delta)$

    Treat $\frac{h^2}{2}f''(\delta)$ as small and ignore the $O(h^2)$ term.
    Then $f(x) + hf'(x) = 0 \implies h = -\frac{f(x)}{f'(x)}$

- Linearize the function $f(x)$ (tangent line) to approximate the curve of $f(x)$ and find the zero of the tangent line.

- Faster than bisection method, but might diverge. Only have **local convergence**.

!!! example
    Q: Use Newton's method to find $1/b > 0$ without using division.

    Solution: Use $f(x) = b-1/x$. Then $1/b$ is the root of $f(x) = 0$.

    $f'(x) = \frac{1}{x^2}$

    $x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)} = x_k - \frac{b - (1/x_k)}{1/x^2} = x_k(2-bx_k)$

    **Extension to matrix inverse**: $X_{k+1} = X_k (2I - AX_k)$

### Error analysis

- Quadratic convergence
- **Assumptions**: $f''(x)$ continuous, $f'(r) \neq 0$

!!! danger
    If the root has multiplicity greater than one, the convergence rate is **linear**.

!!! important "Theorem"
    Let $f''$ be continuous and $r$ satisfy $f(r) = 0, f'(r) \neq 0$. Then ther is **a neighborhood of $r$** and a constant $C$ such that if $x_0$ is in this neighborhood, the successive points generated become steadily closer to $r$ and satisfy 

    $$
    |x_{n+1} - r| \leq C |x_n - r|^2
    $$

??? question "Proof (Taylor's theorem)"
    $e_n = x_n - r$

    Assume: $f''(x)$ continuous, $f'(r) \neq 0$

    $e_{n+1} = x_{n+1} - r = x_n - \frac{f(x_n)}{f'(x_n)} - r = \frac{e_n \cdot f'(x_n) - f(x_n)}{f'(x_n)}$

    By Taylor's theorem,

    $0 = f(r) = f(x_n - e_n) = f(x_n) - e_n f'(x_n) + \frac{1}{2} e_n^2 f''(\delta)$

    $e_n f'(x_n) - f(x_n) = \frac{1}{2} e_n^2 f''(\delta)$

    Substituting,

    $e_{n+1} = \frac{\frac{1}{2} e_n^2 f''(\delta)}{f'(x_n)}$

    As $n\to\infty$, $\delta\to r$ and $x_n\to r$. So $\frac{e_{n+1}}{e_n^2} = \frac{\frac{1}{2} f''(r)}{f'(r)} = C$

### Secant method

Motivation: Calculationg $f'(x)$ is difficult. Why not use an approximation of $f'(x)$ instead?

Approximate $f'(x_n)$ as $\frac{f(x_n)-f(x_{n-1})}{f(x_n)-f(x_{n-1})}$

So the iterative scheme becomes 

$$
x_{n+1} = x_n - \frac{f(x_n)-f(x_{n-1})}{f(x_n)-f(x_{n-1})} f(x_n)
$$

- Slower than Newton's method (proof not required)

!!! warning "Practical considerations for Newton's method"

    - Difficulty in calculating $f'(x)$
    - Failure of the method to converge to the root
    - Overshoot and diverge away from the root
    - Stationary point $f'(x) = 0$
    - Poor initial estimate leading to non-convergence
    - Slow convergence for repeated roots
        - **If the root has multiplicity greater than one, the convergence rate is linear.**
    - Divergent loop

## System of nonlinear equations

$$
\begin{cases} f_1(x_1, x_2) = 0 \\ f_2(x_1, x_2) = 0 \end{cases}
$$

By Taylor's Theorem,

$$
f(x_1^{(k)} + h_1^{(k)}, x_2^{(k)} + h_2^{(k)}) = f(x_1^{(k)}, x_2^{(k)}) + h_1^{(k)} \frac{\partial f}{\partial x_1}(x_1^{(k)}, x_2^{(k)}) +  h_2^{(k)} \frac{\partial f}{\partial x_2}(x_1^{(k)}, x_2^{(k)}) + O(h^2) = 0
$$

Jacobian matrix:

$$
J(X^{(k)}) = \begin{bmatrix}
\frac{\partial f_1}{\partial x_1}(X^{(k)}) & \frac{\partial f_1}{\partial x_2}(X^{(k)})\\
\frac{\partial f_2}{\partial x_1}(X^{(k)}) & \frac{\partial f_2}{\partial x_2}(X^{(k)})
\end{bmatrix}
$$

$$
Jh = -F(X^{(k)})
$$

$$
X^{(k+1)} = X^{(k)} - (J(X^{(k)}))^{-1}F(X^{(k)})
$$

Multi-variable Newton's method still enjoys the property of **quadratic convergence** of the starting point is near the exact solution.

However, a significant weakness is the requirement that, in each iteration:

- a Jacobian matrix (involving $n^2$ partial derivatives) has to be evaluated
- an $n \times n$ linear system involving this matrix must be solved.

### Rescaled simple iteration: Reduce cost of approximating $F'(X^{(k)})^{(-1)}$

- Simplest way: Not to recompute it, just use an initial approximation $F'(X^{(0)})^{(-1)}$ throughout the iteration
- Difficult to guarantee that this condition will be met... Convergence is generally linear
- Or take one Newton step followed by several additional iteration steps without updating the Jacobian matrix.

??? abstract "Generalized secant method (Broyden method) (not required)"
    Begin with a rough guess of the Jacobian, and use successive evaluations of $F$ and its gradient to evaluate the guess of the Jacobian $J$.
    
    Each iteration is far less costly, but more iterations are needed. 

## Fixed point iteration

$$
x_{n+1} = g(x_n)
$$

Finding a zero/root for $f(x)$ $\implies$ Finding a fixed point for $g(x) = x - f(x)$

!!! important "Contractive Mapping Theorem"
    Contractive if there exists $\lambda \in (0, 1)$ s.t. $|F(x) - F(y)| \leq \lambda |x - y|$

    Suppose $F: C \to C$ where $C$ is a closed set, is a contractive mapping. Then $F$ has a _unique_ fixed point.
    
    This fixed point is the limit of every sequence obtained by $x_{n+1} = F(x_n)$ with any starting point $x_0 \in C$

??? question "Proof of convergence of $x_n$"
    Write $x_n = x_0 + (x_1 - x_0) + (x_2 - x_1) + \cdots (x_n - x_{n-1})$

    Need to show that $\sum^{\infty}_{n=1}(x_n - x_{n-1})$ converges. Suffices to show that $\sum^{\infty}_{n=1}|x_n - x_{n-1}|$ converges.

    **Since $F$ is contractive, we have $|x_n - x_{n-1}| = |F(x_{n-1}) - F(x_{n-2})| \leq \lambda |x_{n-1} - x_{n-2}|$**

    This argument can be repeated:

    $|x_n - x_{n-1}| \leq \lambda |x_{n-1} - x_{n-2}| \leq \lambda^2|x_{n-2} - x_{n-3}| \leq \cdots \leq \lambda^{n-1}|x_1 - x_0|$

    Then we have $\sum^{\infty}_{n=1}|x_n - x_{n-1}| \leq \sum^{\infty}_{n=1}\lambda^{n-1}|x_1-x_0| = |x_1-x_0|\sum^{\infty}_{n=1}\lambda^{n-1}=\frac{1}{1-\lambda}|x_1-x_0|$

    So the sum is bounded, hence it converges, and thus the sequence $x_n$ converges for any inital point $x_0$.

??? question "Proof of uniqueness of fixed point"
    Let $x$ and $y$ both be fixed points of $F$.

    Then $|x - y| = |F(x) - F(y)| \leq \lambda |x - y|$

    Since $\lambda < 1$, this forces $|x - y| = 0$. That is, $x = y$.

!!! example
    Q: Prove that the sequence $x_0 = -15, x_{n+1} = 3 - \frac{1}{2}|x_n|$ is convergent.

    A: $|F(x) - F(y)| = |3 - \frac{1}{2}|x| - 3 + \frac{1}{2}|y|| = \frac{1}{2}||y|-|x|| \overset{\mathrm{Triangle~inequality}}{\leq} \frac{1}{2}|y-x|$

!!! important "Useful convergence criteria"
    If $F \in C[a, b]$ such that $a \leq F(x) \leq b$ for all $x \in [a, b]$, then $F$ has a fixed point in $[a, b]$. ($F$ maps a closed set onto itself)
    
    In addition, if $F'(x)$ exists in $(a, b)$ and there exists a positive constant $M < 1$ such that $|F'(x)| \leq M < 1$ for all $x \in (a, b)$, the fixed point is unique. ($F$ is contractive)

    Smaller slope will lead to faster convergence!

??? question "Proof"
    Existence of fixed point: Intermediate value theorem on $F(a) - a > 0$ and $F(b) - b < 0$.

    Unique fixed point: Mean value theorem on two fixed points, must exist slope $1$ between them, but this contradicts the slope bound condition.

### Error analysis

Order of convergence is $m$, where $F^{(k)}(x^*) = 0$ for $1 \leq k \leq m-1$ but $F^{(k)}(x^*) \neq 0$

??? question "Proof"
    $e_{n+1} = x_{n+1} - x^* = F(x_n) - F(x*) = F'(x^*)(x_n - x^*) + \cdots + \frac{F^{(m-1)}(x^*)}{(m-1)!}(x_n - x^*)^{m-1} + \frac{F^{(m)}(\delta_n)}{m!}(x_n - x^*)^{m} = \frac{F^{(m)}(\delta_n)}{m!}e_n^m$

    So $\lim_{n\to\infty}\frac{|e_{n+1}|}{|e_n|^m} = \frac{1}{m!}|F^{(m)}(x^*)|$