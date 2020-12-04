# Chapter 5: Estimation

## Rao-Blackwell Theorem

Define $\rho^* = \rho^*(T) = \mathbb{E}[\rho(X) | T]$

1. $\rho^*$ is an estimator of $\psi(\theta)$ (it is a function of data, not depending on $\theta$)
2. $\mathbb{E}[L(\theta, \rho^*)] \leq \mathbb{E}[L(\theta, \rho)]$ (MSE)
3. If $\rho$ unbiased and $T$ complete for $\theta$, then
    - $\rho^* = \rho^*(T)$ is the **only** function of $T$ which is unbiased for $\psi(\theta)$
    - $\forall$ **unbiased** estimator $S$ of $\psi(\theta)$, $\mathbb{E}[L(\theta, \rho^*)] \leq \mathbb{E}[L(\theta, S)] \forall \theta$ (UMVU estimator under squared loss)

!!! success "Important implications"
    Provide constructive way to find "optimal" (smallest risk & unbiased) estimator (e.g. UMVU estimator) if complete sufficient statistic $T$ exists. Either:

    - Choose any unbiased estimator $\rho$, modify it to $\rho^* = \mathbb[\rho|T]$
    - Find function $g(T)$ which is unbiased, then just use $g(T)$

    They will get the same result from 3.2.

    Also, Rao-Blackwell highlights an important use of complete sufficient statistic. This alone justifies our need to study the concept of **completeness**.

!!! quote "Exam example: 2017"
    ... some calculations of MLE and complete sufficient statistic

    The mle $\hat{\theta}$ is a function of $X-Y$ according to (previous part). If it were unbiased for $\theta$, then it must be identical to the UMVU estimator by Rao-Blackwell Theorem. That the two estimates differ on the data, which has a positive probability to be observed, implies that $\hat{\theta}$ cannot be unbiased for $\theta$.


!!! example
    ![ch5 rao blackwell eg](../images/ch5_rao_blackwell_eg.png)

## Fisher information

!!! important "Motivation"
    Compare experiments and their likelihood functions. It is better to have a likelihood function that is "sharper" at the true value.

    Sharper means the negative of the second derivative of loglikelihood!

!!! important "Theorem: Information matrix"
    ![ch5 information matrix](../images/ch5_information_matrix.png)

!!! success "Tip"
    ![ch5 information matrix tip](../images/ch5_information_matrix_tip.png)

## Information inequality, CRLB

!!! warning "Regularity conditions"
    - Sample space $S$ must not depend on $\theta$
        - $U[0, \theta]$ ruled out
    - Parameter space $\theta$ contains open rectangle
        - Discrete $\theta$ ruled out
    - $f(x|\theta)$ twice differentiable w.r.t. $\theta$
        - $\text{Binomial}(\theta, p)$ ruled out, $\text{Binomial}(n, \theta)$ ok

??? question "What if the regularity assumptions are violated?"
    ![ch5 violate regularity](../images/ch5_violate_regularity.png)

!!! note "CRLB"
    CRLB: $\text{Var}_\theta(\hat{\theta}_U) \geq \frac{1}{I(\theta)}$

!!! important "Theorem"
    ![ch5 information inequality](../images/ch5_information_inequality.png)

## Maximum likelihood estimator (MLE)

Obtained by solving likelihood equations: $U(\theta) = \frac{\partial S_X(\theta)}{\partial \theta} = 0$

- MLE is function of minimal sufficient statistic $T$, but $\hat{\theta}$ may not be sufficient.
- MLE may not be unbiased
- $\hat{\theta}$ mle of $\theta$, then $\psi(\hat{\theta})$ mle of $\psi(\theta)$ 

## Large sample properties of MLE

- $\hat{\theta}_n \to \theta_0$ in probability
- $\sqrt{n} (\hat{\theta}_n - \theta_0) \to N(0, \mathscr{I}(\theta_0)^{-1})$ **(very important!!!)**
    - $\mathscr{I}(\theta_0) = \lim_{n\to\infty} n^{-1}I(\theta_0)$
    - $\hat{\theta}_n \sim N(\theta_0, I(\theta_0)^{-1})$ (inverse of Fisher information matrix)
- $n^{-1/2} U(\theta_0) \to N(0, \mathscr{I}(\theta_0))$

!!! important "Asymptotic distribution of $\hat{\theta}_n$"
    ![ch5 asymptotic distribution](../images/ch5_asymptotic_distribution.png)

!!! warning
    Would **not** help if $\psi(\hat{\theta}_n)$ is changed to $\mathbb{E}[\psi(\hat{\theta}_n) | \text{sufficient statistic}]$, as already $\psi(\hat{\theta}_n) = \mathbb{E}[\psi(\hat{\theta}_n) | \text{sufficient statistic}]$

!!! note
    In practice, to obtain the asymptotic distribution from data, we can calculate $\hat{\theta}_n$ and $I(\hat{\theta}_n)^{-1}$ and treat them as $\theta_0$ and $I(\theta_0)^{-1}$ respectively.
