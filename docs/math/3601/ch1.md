# Chapter 1: Preliminaries

## Errors

Absolute error: $|x - \tilde{x}|$

Relative error: $\frac{|x - \tilde{x}|}{|x|}$

## Important theorems

!!! important "Rolle's Theorem"

!!! important "Intermediate Value Theorem"

!!! important "Mean-Value Theorem"

!!! important "Taylor's Theorem"
    $R_n$: Truncation error estimate

    2D Case

!!! important "Triangle inequality"
    $||a|-|b|| \leq |a-b|$

!!! important "AM GM inequality"
    $\frac{x+y}{2} \geq \sqrt{xy}$

## Condition number


$$
CN = \frac{\frac{f(x) - f(x_0)}{f(x_0)}}{\frac{x-x_0}{x_0}} \approx \frac{x_0 f'(x_0)}{f(x_0)}
$$

## Order of Convergence

$$
\lim_{n\to\infty} \frac{|x_{n+1}-x^*|}{|x_n-x^*|^p} = \beta
$$

Use L'Hopital rule?

## Big O

$O(h^2)$