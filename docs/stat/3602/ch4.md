# Chapter 4: Sufficiency and Likelihood

!!! question "Fundamental problem in parametric statistical inference: How to infer about $\theta$ from evidence $X$, $X \sim f(x|\theta)$?"

!!! question "Can observation $X$ help at all?"

**Statistic**: a "partition" of sample space and "carrier" of information

## Sufficient statistic

!!! important "Definition: $T = T(X)$ sufficient for $\theta$ if distribution of $X$ conditional on $T$ does not depend on $\theta$"

    i.e. given knowledge of $T(X)$, no extra information relevant to $\theta$ can be gained by knowing further details of $X$

    info contained in $T(X)$ about $\theta$ = info contained in $X$ about $\theta$

    Sufficiency principle: Under same experimental setting, if $T(X) = T(X')$, then inference from $X_A$ $\equiv$ Inference from $X_B$.

    Outcomes in same cell of partition (same $T(X)$) share proportional likelihood function $l(\theta)$.

Sufficient statistic is not unique.

## Likelihood function

$l_X(\theta) \propto f(X|\theta)$

- Which $\theta$ is more likely
- Proportionality constant to scale the sum to $1$.
- "posterior probability function" derived from non-informative prioir. but the area of $l_X(\theta)$ may be $\infty$

### Loglikelihood function

$S_X(\theta) = \ln l_x(\theta)$

$S_X(\theta) = \ln l_X(\theta) = \ln \prod^n_{i=1} l_{X_i}(\theta) = \sum^n_{i=1} S_{X_i}(\theta)$

- sum of independent terms, good for **Central Limit Theorem**, **WLNN**, **SLNN**, ...

!!! important "Theorem"
    ![sufficiency theorem](../images/ch4_sufficiency_thm.png)

    ??? question "Proof (proportional likelihood function)"
        ![prop likelihood](../images/ch4_prop_likelihood.png)

    ??? question "Proof (factorization theorem)"
        ![factorization](../images/ch4_factorization.png)

    ??? question "Proof (factorization is sufficient)"
        ![factorization is sufficient](../images/ch4_factorization_is_sufficient.png)
    
??? example
    ![ch4 eg sufficiency](../images/ch4_eg_sufficiency.png)

## Minimal sufficient

$T$ is sufficient for $\theta$ and is function of any sufficient statistic for $\theta$

$\forall X, X', T(X) = T(X') \Leftrightarrow l_X(\theta) \propto l_{X'}(\theta)$

!!! warning "Remark"
    $T(X)$ sufficient for $\theta$:

    iff $\forall X, X', T(X) = T(X') \Rightarrow l_x(\theta) \propto l_{x'}(\theta)$

    One direction only!!

??? example "Example (why $\max{X_i}$ is minimal sufficient for $U[0, \theta]$)"
    ![ch4 minimal sufficient eg](../images/ch4_minimal_sufficient_eg.png)

## Sufficiency of exponential family

??? important "Theorem: Natural statistic $T$ is sufficient."
    ![ch4 sufficient exponential](../images/ch4_sufficient_exponential.png)

??? important "Theorem: If $\Pi$ is not contained in any affine hyperplane in $\mathbb{R}^k$, then $T(X)$ is minimal sufficient for $\theta$ (NOT single point in $\mathbb{R}_1$, NOT straight line in $\mathbb{R}_2$"
    ![ch4 affine hyperplane](../images/ch4_affine_hyperplane.png)

??? example "Example: 2D rectangle minimal sufficient"
    ![ch4 affine eg1](../images/ch4_affine_eg1.png)

??? example "Example: 2D straight line not minimal sufficient"
    ![ch4 affine eg2](../images/ch4_affine_eg2.png)

??? example "Example: 2D curve minimal sufficient"
    ![ch4 affine eg3](../images/ch4_affine_eg3.png)

??? important "Theorem: If $\Pi$ contains open rectangle in $\mathbb{R}^k$, then natural statistic $T(X)$ is complete sufficient for $\pi$."
    Open rectangle: Open interval in 1d, Open rectangle in 2d

    !!! question "Proof"
        ![ch4 open rectangle proof](../images/ch4_open_rectangle_proof.png)

!!! quote
    The natural parameter space for $\ln{\theta}$ is $(-\infty, \infty)$, which contains an open interval. Thus $T$ is a complete sufficient statistic for $\theta$.

??? example "2D rectangle complete"
    ![ch4 open rectangle eg1](../images/ch4_open_rectangle_eg1.png)

??? example "2D curve not complete"
    ![ch4 open rectangle eg2](../images/ch4_open_rectangle_eg2.png)

!!! warning
    If not contain open rectangle, it is possible that it is still complete sufficient

    ??? example
        ![ch4 open rectangle eg3](../images/ch4_open_rectangle_eg3.png)

!!! warning
    Complete statistic must be minimal sufficient
    
    But minimal sufficient statistic may not be complete

    Not complete statistic can still be minimal sufficient

## Likelihood principle

If $l_A(\theta) \propto l_B(\theta)$, then $\text{Inference from A } \equiv \text{ Inference from B}$ 

## Complete sufficient

Definition: $T(X)$ sufficient, $\mathbb{E}_\theta[g(T)] = 0 \forall \theta \Rightarrow \mathbb{P}(g(T) = 0) = 1 \forall \theta$

!!! important "Theorem (Lehmann-Scheffe): $T(X)$ complete sufficient for $\theta$ $\Rightarrow$ $T(X)$ minimal sufficient for $\theta$"
    ??? question "Proof"
        ![ch4 lehmann scheffe proof](../images/ch4_lehmann_scheffe_proof.png)

!!! example
    ![ch4 complete eg1](../images/ch4_complete_eg1.png)
    
!!! example
    ![ch4 complete eg2](../images/ch4_complete_eg2.png)

!!! example
    ![ch4 complete eg3](../images/ch4_complete_eg3.png)

## Remark on proportionality

!!! warning
    ![ch4 propto remark](../images/ch4_propto_remark.png)