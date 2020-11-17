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
    
!!! example
    ![ch4 eg sufficiency](../images/ch4_eg_sufficiency.png)


