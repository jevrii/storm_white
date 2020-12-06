# Chapter 6: Hypothesis testing

## Some definitions

!!! abstract
    Test function: $\phi(X) = \mathbb{P}(\text{reject }H_0 | X)$

    - If we observe $X$, then reject $H_0$ with probability $\phi(X)$

    Power function: $\omega(\theta) = \mathbb{E}_\theta [\phi(X)]$

    - Re: Type I error, 1 - Type II error
    - Probablity of reject given the true $\theta$

    Size: $\sup_{\theta \in \Theta_0} \omega(\theta) = \alpha$

    - Probability of rejecting true $H_0$

    Power: $\omega(\theta)$ at $\theta \in \Theta_1$

    p-value: $\sup_{\theta \in \Theta_0} \mathbb{P}_\theta (T(X) \geq T(x))$

    - $x$ is observation, $T(X)$ is random and has a known distribution
    - Smallest size of LR test that rejects $H_0$ (refer to ranking example)

    Test is **unbiased** for size $\alpha$ if $\sup_{\theta \in \Theta_0} \omega(\theta) = \alpha$, and $\omega(\theta) \geq \alpha \forall \theta \in \Theta_1 \setminus \Theta_0$ 

    UMP: Size $\leq \alpha$, and power is the greatest at $\Theta_1$ for all tests. If it exists, it is unbiased.

    UMPU: Size $\leq \alpha$, and power is the greatest at $\Theta_1$ for all unbiased tests

## Note on critical region and p-value

p-value: $\mathbb{P}_{\theta_0} (T(X) \geq T(x))$

Critical region: $\mathbb{P}_{\theta_0} (T(X) > c) = \alpha$

## Likelihood ratio test

![ch6 likelihood ratio test](../images/ch6_likelihood_ratio_test.png)

!!! important "Neyman Pearson Lemma"
    For simple hypothesis, likelihood ratio test is most powerful among all tests with size $\alpha$

!!! example
    ![ch6 lrtest eg1](../images/ch6_lrtest_eg1.png)

## UMP test under monotone likelihood ratio (mlr property)

See if the likelihood ratio is monotonic respect to some statistic $T(X)$

Need the likelihood ratio to be increasing. So might need to take $-T(X)$.

!!! important "Theorem"
    Test function $T(X) > t_0$ is LR test, unbiased, UMP among tests of size $\leq \sup_{\theta \in \Theta_0} \mathbb{E}_\theta [\phi_0(X)]$

!!! note "Note for exponential family: The sign of $T$ follows sign of test."
    ![ch6 mlr eg1](../images/ch6_mlr_eg1.png)

??? example
    ![ch6 mlr eg2](../images/ch6_mlr_eg2.png)

## Conditional test for exponential family

Only want to know single parameter and ignoring others

!!! important "Theorem"
    ![ch6 conditional 1](../images/ch6_conditional_1.png)

??? example "Example: Normal"
    Refer to slides

    fuck

??? example "Example: Poisson"
    ![ch6 conditional eg2](../images/ch6_conditional_eg2.png)

## UMPU test for exponential family

Two sided: $H_0: \theta \in [\theta_1, \theta_2]$ vs $H_1: \theta \notin [\theta_1, \theta_2]$

Actually can just refer to the conditional case

??? abstract "Summary"
    ![ch6 test summary](../images/ch6_test_summary.png)

??? imporant "Theorem (UMPU interval two-sided)"
    ![ch6 umpu exponential 1](../images/ch6_umpu_exponential_1.png)

??? imporant "Theorem (UMPU Single point two sided: $H_0: \theta = \theta_0$ vs $H_1: \theta \neq \theta_0$)"
    ![ch6 umpu exponential 2](../images/ch6_umpu_exponential_2.png)

??? example "Example 1: Interval"
    ![ch6_umpu_eg1_1.png](../images/ch6_umpu_eg1_1.png)

    ![ch6_umpu_eg1_2.png](../images/ch6_umpu_eg1_2.png)

??? example "Example 2: Point"
    ![ch6_umpu_eg2.png](../images/ch6_umpu_eg2.png)

## Confidence set

Invert the hypothesis test: Find the $\theta_0$ such that $H_0$ will not be rejected.

Example: Test $\theta = k$ vs $\theta > k$. 

Then during hypothesis testing, we will obtain the p-value $\mathbb{P}(T_{random} > t_{obs} | \theta = k)$.

So to obtain the confidence bound, we solve $\mathbb{P}(T_{random} > t_{obs} | \theta = k) \geq \alpha$ on $k$.

!!! note
    ![ch6 confidence set 1](../images/ch6_confidence_set_1.png)

??? example
    ![ch6 confidence set eg1](../images/ch6_confidence_set_eg1.png)

??? quote "2018 exam"
    ![ch6_2018_exam_eg1.png](../images/ch6_2018_exam_eg1.png)

## Generalized Likelihood Ratio (GLR) Test

!!! Warning
    As this is a test based on large-sample theory, the size may be greater than $\alpha$ if $n$ is too small...

1. Find the degree of freedom $s$ for $H_1$, $r$ for $H_0$. $H_0$ is more constrained than $H_1$ (e.g. $\mu = \lambda$)
2. $G = 2(\ln(\ell_X(\hat{\theta}) - \ln(\ell_X(\tilde{\theta})))) > \chi^2_{r-s}(\alpha)$

Samples are independent, not necessarily identically distributed
Samples are independent, not necessarily identically distributed

!!! important "Theorem"
    ![ch6 glr v2](../images/ch6_glr_v2.png)

??? example
    ![ch6 glr 2](../images/ch6_glr_2.png)

!!! example
    ??? note "Part 1: Find MLE and constrained MLE"
        ![ch6 glr eg1 1](../images/ch6_glr_eg1_1.png)

    ??? note "Part 2: Find $G = 2(\ln(\ell_X(\hat{\theta}) - \ln(\ell_X(\tilde{\theta})))) > \chi^2_{r-s}(\alpha)$"
        ![ch6 glr eg1 2](../images/ch6_glr_eg1_2.png)
