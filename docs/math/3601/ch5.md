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

### Trapezoidal Rule

$$
I(f) \approx \frac{b-a}{2}(f(a) + f(b))
$$

Exact for linear function.

### Simpson's rule

Quadratic Newton-Cotes formula.

$$
\int^b_a f(t) dt = \frac{b-a}{6}(f(a) + 4f((a+b)/2)+f(b))
$$

Exact for **cubic** function. (some magic here??)

!!! danger "The order of accuracy is not continuous!"
    Re: Even nodal point, Odd nodal point

### Method of undetermined coefficients

Solve 

$$
\int^b_a f(x) dx = \sum^n_{j=0} A_j f(x_j)
$$

Take $f(x) = 1, f(x) = x, f(x) = x^2$. Then solve the system of linear equations.

## Composite rules

### Composite Simpson's rule

$$
\frac{b-a}{3n}(f(x_0) + 4f(x_1) + f(x_2) + f(x_2) + 4f(x_3) + f(x_4) + f(x_4) + \cdots + 4f(x_{n-1}) + f(x_n))
$$