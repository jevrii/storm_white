# Chapter 3: Exponential Families

## Exponential family

!!! important "Theorem: Exponential family"

    $$
    f(x|\theta) \propto h(x) \exp(\sum^k_{j=1} \pi_j(\theta) t_j(x)), x \in S \text{ (free of } \theta \text{)}
    $$

    $$
    f(x|\pi) \propto h(x) \exp(\sum^k_{j=1} \pi_j t_j(x)), x \in S \text{ (free of } \theta \text{)}
    $$

    $$
    f(x|\pi) = C(\pi) h(x) \exp(\sum^k_{j=1} \pi_j t_j(x)) \text{ (normalizing constant)}
    $$

    Note: Uniform distribution is not exponential family form, because $x \in [0, \theta]$, the sample space depends on $\theta$.

!!! example
    ![exponential form of common pdf](../images/ch3_exponential_form_pdf.png)

## Natural parameter space

- Natural parameter space: $\Pi = \{ \pi = (\pi_1(\theta), \pi_2(\theta), \cdots, \pi_k(\theta)) : \theta \in \Theta \} \subset \mathbb{R}^k$
- Natural parameter $\pi = (\pi_1, \cdots, \pi_k) \in \Pi$
- Natural statistic $(t_1(X), \cdots, t_2(X))$ for $X \in S$ drawn from $f(\cdot | \theta)$
    - e.g. $T = \sum^n_{i=1}X_i$, $T = (\sum^n_{i=1}X_i^2, \sum^n_{i=1}X_i)$

Natural parameter space $\Pi$ can be different from parameter space $\Theta$

!!! example
    $\text{Poisson}(\lambda)$: $\theta = \lambda \in \Omega = (0, \infty), \pi_1(\lambda) = ln(\lambda) \in \Pi = (-\infty, \infty)$

    $\text{N}(\mu, \sigma^2)$: $\begin{align*}\theta = (\mu, \sigma^2) \in \Omega = (-\infty, \infty) \times (0, \infty),\\ (\pi_1(\lambda), \pi_2(\lambda)) = \left(-\frac{1}{2\sigma^2}, \frac{\mu}{\sigma^2}\right) \in \Pi = (-\infty, 0) \times (-\infty, \infty)\end{align*}$

    $\text{N}(|\theta|, |\theta|+\theta^2)$: $\theta \in \Omega = (-\infty, -1] \cap [1, \infty)$, $(\pi_1(\theta), \pi_2(\theta)) = \left(\frac{1}{2(|\theta + \theta^2)}, \frac{1}{1+|\theta|}\right) \in \Pi = \{\pi_1 = \frac{\pi_2^2}{2(\pi_2-1)}\}$

![exponntial sample iid](../images/ch3_exponential_sample_iid.png)

!!! important "Theorem: Distribution of natural statistic $T = (t_1(X), \cdots, t_k(X))$ has exponential family form with natural parameter also given by $\pi \in \Pi$"
    ![ch3 proof](../images/ch3_proof.png)

??? example
    ![ch3 poisson eg](../images/ch3_poisson_eg.png)

## Conditional distribution of natural statistic

!!! theorem "Conditional distribution of natural statistic"
    Distribution of $T_1$ conditional on $T_2$ has exponential family form with natural parameter $\pi_1$

    Can get rid of $\pi_2$

!!! question "Proof"
    $\text{Conditional} = \frac{\text{Joint of } T_1, T_2}{\text{Marginal of } T_2}$

    $g(t_1 | t_2, \pi) = \frac{g(t_1, t_2 | \pi)}{g(t_2 | \pi)} \propto g(t_1, t_2 | \pi) \propto h^* (t_1, t_2) exp(\pi_1 t_1 + \pi_2 t_2) \propto h^{**} exp (\pi_1 t_1)$

So we can try to formulate $\pi_1$ such that it is useful to us.

!!! example
    ![ch3 conditional eg](../images/ch3_conditional_eg.png)

