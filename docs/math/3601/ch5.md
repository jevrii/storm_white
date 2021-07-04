# Chapter 5: Numerical Integration

!!! quote
    The differential calculus is a **science**; the integral calculus is an **art**.

$$
I(f) \approx Q(f) := \sum^n_{j=0} A_j f(x_j)
$$

Q: How to determine weights $A_j$? How to measure the error?

## Integration by interpolation polynomials

$p(x) = \sum^n_{i=0} f(x_i)l_i(x)$

$I(f) \approx \sum^n_{i=0}f(x_i) \int^b_a l_i(x) dx$

When the $x_i$'s are equally spaced, then this formula is called the **Newton-Cotes** formula of order $n$.

!!! important "Theorem"
    The order of accuracy is $n$ if it is exact for all $f \in \Pi_n$ but there will be a nonzero error term for some degree $(n+1)$ polynomial.

    The formula is exact for all $f \in \Pi_n$ if and only if it is derived from the integration of the polynomial (with minimal degree) interpolating $f(x_0), f(x_1), \cdots, f(x_n)$.

    $A_i = \int^b_a \ell_i(x) dx$

### Trapezoidal Rule

$$
I(f) \approx \frac{b-a}{2}(f(a) + f(b))
$$

Composite trapezoidal rule:

$$
\begin{align*}
\frac{h}{2}\left(f(a)+f(b)+2\sum^{n-1}_{i=1}f(x_i)\right)\\
x_i = a+ih, h=(b-a)/n
\end{align*}
$$

Exact for linear function. Order of accuracy is one.

!!! note "Error analysis of trapezoidal rule"
    If $f \in C^2$, then the error of the trapezoidal's rule is $-\frac{1}{12}f''(\xi)(b-a)^3$, where $\xi \in (a, b)$

!!! note "Error analysis of composite trapezoidal rule"
    If $f \in C^2$, then the error is $-\frac{h^3}{12}\sum^{n-1}_{i=0}f''(c_i) = -\frac{h^2}{12}(b-a)f''(c)$, where $c_i \in (x_i, x_{i+1})$ and $c \in [a, b]$.


### Simpson's rule

Quadratic Newton-Cotes formula.

$$
\int^b_a f(t) dt = \frac{b-a}{6}(f(a) + 4f((a+b)/2)+f(b))
$$

Composite Simpson's rule:

$$
\frac{b-a}{3n}(f(x_0) + 4f(x_1) + f(x_2) + f(x_2) + 4f(x_3) + f(x_4) + f(x_4) + \cdots + 4f(x_{n-1}) + f(x_n))
$$

!!! note "Error analysis of Simpson's rule"
    if $f \in C^4$, then the error of the Simpson's rule is $-\frac{1}{90}f^{(4)}(\xi)h^5$, where $h=(b-a)/2$ and $\xi \in (a, b)$

Exact for **cubic** function (some magic here??). Order of accuracy is three.

!!! danger "The order of accuracy is not continuous!"
    Re: Even nodal point, Odd nodal point

!!! note "Linearity of integration and interpolation"
    For Simpson's rule, we can verify the cubic exact by putting $f(x) = x^3$. **By the linearity of integration and the interpolation**, the formula is exact for $p \in \Pi^3$.
    
    It is not exact when $f(x) = x^4$

### Method of undetermined coefficients

Solve ($x_j, a, b, n$ given)

$$
\int^b_a f(x) dx = \sum^n_{j=0} A_j f(x_j)
$$

Take $f(x) = 1, f(x) = x, f(x) = x^2, \cdots, x_n$. Then solve the system of linear equations. It is equivalent to taking a basis in the linear space of functions, so the functions can be "a linear combination of nonlinear functions" (e.g. $f(x) = ae^x + b \cos(\pi x/2)$)

## Gaussian Quadrature 

!!! success "Motivation: Choose both the weights $A_i$ and nodes $x_i$ freely, then we can achieve exact for $f \in \Pi_{2n+1}$"

1. Find the q-orthogonal polynomial and the roots
2. Solve the weigths

!!! example
    ![ch5_gaussian_quadrature_eg1.png](../images/ch5_gaussian_quadrature_eg1.png)

!!! important "Gaussian quadrature"
    ![ch5_gaussian_quadrature.png](../images/ch5_gaussian_quadrature.png)

??? abstract "Gram-Schmidt Process"
    ![ch5_gram_schmidt](../images/ch5_gram_schmidt.png)

??? question "Proof"
    !!! important "Definition: w-orthogonal"
        ![ch5_orthogonal_def.png](../images/ch5_orthogonal_def.png)

    ![ch5_gaussian_quadrature_proof.png](../images/ch5_gaussian_quadrature_proof.png)

