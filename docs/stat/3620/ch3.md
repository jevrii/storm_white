# Chapter 3: k-Sample methods

When there are more than two treatments:

- Do an overall comparison to determine whether or not there are differences among the $k$ samples
- If the differences are found, we perform multiple comparisons to further detect which treatments differ from each other.

Test between

- $H_0: F_1(x) = F_2(x) = \cdots = F_k(x)$ for all $x$, and
- $H_1: F_i(x) \leq F_j(x)$ or $F_i(x) \geq F_j(x)$ for at least one pair with strict inequality holding for at least one $x$.

Alternatively, let $F_i(x) = F(x - \mu_i)$. Then the hypothesis become

- $H_0: \mu_1(x) = \mu_2(x) = \cdots = \mu_k(x)$ for all $x$, and
- $H_1: \mu_i$ are not all equal


## Parametric method based on Normality Assumption

$X_{ij} = \mu_{i} + \epsilon_{ij}$

$SST = SSB + SSE$

F-test: $F = \frac{SSB/(k-1)}{SSE/(N-k)} \sim F(k-1, N-k)$ under $H_0$, large $F$ values favor the rejection of $H_0$.

However, if we are unwilling to assume that the population distributions are normal, we may carry out a permutation version of F-test instead.

## Permutation F-test

1. Compute the F-statistic for the observed data, denoted $F_{obs}$.
2. For each of the $\frac{N!}{n_1!n_2! \cdots n_k!}$ permutations (or a random sample of the permutations), compute the F-statistic
3. Calculate the p-value as the fraction of the F's that are greater than or equal to $F_{obs}$.

Can also perform permutation test based on $SSB = \sum^k_{i=1}n_i(\bar{X_i} - \bar{X})^2$ or $SSX = \sum^k_{i=1} n_i \bar{X_i}^2$ as the test statistic.

!!! note
    Note that **for large samples, the permutation
    distribution may be approximated by the F-distribution. In other words, the ANOVA
    test may well approximate the permutation F-test.**

    This is easy to conprehend since sample
    means will be normally distributed for large samples, according to the Central Limit Theorem.

!!! warning
    In practice, the ANOVA test is relatively robust to departures from the normality assumption
    but it is sensitive to departures from the **constant variance assumption**. [link](http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r)

??? example "Example 3.1"
    Three preservatives were compared to see their ability to inhibit the growth of
    bacteria. Samples were treated with one of the preservatives at the same time or left untreated
    as a control. After 48 hours, bacteria counts were made. The logarithms of the counts (Why?)
    were computed to meet ANOVA assumptions, recorded as the table:

    ```r
    ctrl = c(4.302, 4.017, 4.049, 4.176)
    pre1 = c(2.021, 3.190, 3.250, 3.276, 3.292, 3.267)
    pre2 = c(3.397, 3.552, 3.630, 3.578, 3.612)
    pre3 = c(2.699, 2.929, 2.785, 2.176, 2.845, 2.913)
    x = c(ctrl, pre1, pre2, pre3)
    grps = rep(c("ctrl","pre1","pre2","pre3"), c(4,6,5,6))
    # grps = rep(1:4, c(4,6,5,6))

    dat <- data.frame(x, grps)
    aggregate(x ~ grps, dat, mean)    # compute group means

    # one-way ANOVA model
    onewayanova <- (aov(x ~ factor(grps)))
    summary(onewayanova)

    str(summary(onewayanova))
    Fobs <- summary(onewayanova)[[1]][1,4]
    Fobs

    # Permutation-F test
    R <- 9999
    reps <- numeric(R)

    for (i in 1:R) {
    permdat <- data.frame(x=dat$x,
    grps=sample(dat$grps, size=dim(dat)[1], replace=FALSE))
    permanova <- aov(x ~ factor(grps), data=permdat)
    reps[i] <- summary(permanova)[[1]][1,4] # extract the ith permuted F statistic
    }
    pvalue <- mean(c(Fobs, reps) >= Fobs) # always an upper-tailed test
    pvalue

    hist(c(Fobs,reps), main="Null distribution of F statistic")
    abline(v=Fobs, lty=2 ,col="red")
    ```

    ```
    > # one-way ANOVA model
    > onewayanova <- (aov(x ~ factor(grps)))
    > summary(onewayanova)
    Df Sum Sq Mean Sq F value Pr(>F)
    factor(grps) 3 5.476 1.8254 17.66 1.82e-05 ***
    Residuals 17 1.757 0.1034
    –-
    Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    >
    > str(summary(onewayanova))
    List of 1
    $ :Classes ‘anova’ and ’data.frame’: 2 obs. of 5 variables:
    ..$ Df : num [1:2] 3 17
    ..$ Sum Sq : num [1:2] 5.48 1.76
    ..$ Mean Sq: num [1:2] 1.825 0.103
    ..$ F value: num [1:2] 17.7 NA
    ..$ Pr(>F) : num [1:2] 1.82e-05 NA
    - attr(*, "class")= chr [1:2] "summary.aov" "listof"
    > Fobs <- summary(onewayanova)[[1]][1,4]
    > Fobs
    [1] 17.65729
    ```

    One-way ANOVA has F-statistic = 17:66 with p-value = 0:000018 while the permutation F-test
    based on R = 10000 random permutations gives the exact p-value = 0:0001 (is this OK?).

    Close agreement is found between the two p-values. Note that **for large samples, the permutation
    distribution may be approximated by the F-distribution. In other words, the ANOVA
    test may well approximate the permutation F-test.** This is easy to conprehend since **sample
    means will be normally distributed for large samples, according to the Central Limit Theorem.**

## Kruskal-Wallis Test (ANOVA based on ranks)

!!! important "Objective"
    Obtaining a statistic based on **ranks**, similar to Wilcoxcon test.

    Approximating permutation test using chi-square, useful for large sample size.

Rank the data combining all treatments. Let $\bar{R_i}$ denote the mean of the ranks of treatment $i$'s observations, assuming the data have no ties.

The **Kruskal-Wallis statistic** is defined as:

$$
KW = \frac{12}{N(N+1)}\sum^k_{i=1}n_i\left(\bar{R_i}-\frac{N+1}{2}\right)^2 = \frac{12}{N(N+1)} \sum^k_{i=1} n_i (\bar{R_i})^2 - 3(N+1)
$$

- where $\frac{N+1}{2}$ is the overall mean of all ranks and $\sum^k_{i=1}n_i\left(\bar{R_i}-\frac{N+1}{2}\right)^2$ can be viewed as the $SSB$ on ranks. (recall $SSB = \sum^k_{i=1}n_i(\bar{X_i} - \bar{X})^2$)
- $\frac{12}{N(N+1)}$ is a constant scaling factor to make the $KW$ statistic to follow **approximately** a chi-square distribution under $H_0$, i.e. $KW \sim \chi^2(k-1)$. (note that $E[\sum^k_{i=1}n_i\left(\bar{R_i}-\frac{N+1}{2}\right)^2] = \frac{N(N+1)}{12}(k-1)$, we want $E[KW] = k-1$.
    - $E[\left(\bar{R_i}-\frac{N+1}{2}\right)^2] = Var(\bar{R_i}) + (E[\bar{R_i}] - \frac{N+1}{2})^2 = \frac{1}{n_i^2}\frac{n_i(N-n_i)\sigma^2}{N-1}$
    - Can get the result from p.20 Expectation and Variance of Wilcoxon statistic. But here we are using average rank instead of rank sum, so divide by ${n_i^2}$
    - $\sigma = \frac{(N-1)(N+1)}{12}$

The "exact" p-value of the $KW$ test can be obtained by using the permutation method.

Reject when the upper $\alpha$% critical value from the chi-square table is less than $KW$.

!!! note
    The permutation $KW$ test tends to have smaller critical values than the chi-square critical values. Hence, if we use the chi-square approximation, rejection of the null hypothesis will assure the same result with the exact $KW$ test, and the exact p-value will be smaller.

??? example "Example 3.2"
    Continue to use the data from example 3.1, but this time the sample size $N = 21$ is moderately large.

    $KW = \cdots = 17.14$, $KW \sim \chi^2(3)$ under $H_0$, the upper 1% critical value from the chi-square table is $11.34$. Thus the test is significant at the 1% level and $H_0$ should be rejected.

    ```r
    ctrl = c(4.302, 4.017, 4.049, 4.176)
    pre1 = c(2.021, 3.190, 3.250, 3.276, 3.292, 3.267)
    pre2 = c(3.397, 3.552, 3.630, 3.578, 3.612)
    pre3 = c(2.699, 2.929, 2.785, 2.176, 2.845, 2.913)
    x = c(ctrl, pre1, pre2, pre3)
    grps = rep(c("ctrl","pre1","pre2","pre3"), c(4,6,5,6))
    # grps = rep(1:4, c(4,6,5,6))

    dat <- data.frame(x, grps)
    aggregate(x ~ grps, dat, mean)    # compute group means

    # K-W test 
    kruskal.test(x ~ grps, data=dat)

    # Alternative way of computing the K-W test statistic
    rank.x = rank(x)
    summary(aov(rank.x ~ factor(grps)))
    SSB = summary(aov(rank.x ~ factor(grps)))[[1]][1,2]
    SR2 = var(rank.x)
    KW = SSB/SR2

    # Exact p-value using permutation test
    R <- 9999
    reps <- numeric(R)

    for (i in 1:R) {
    permdat <- data.frame(x=rank.x,
    grps=sample(dat$grps, size=dim(dat)[1], replace=FALSE))
    permkw <- aov(x ~ factor(grps), data=permdat)
    reps[i] <- summary(permkw)[[1]][1,2] / SR2
    }
    pvalue <- mean(c(KW, reps) >= KW) # always an upper-tailed test
    pvalue

    hist(c(KW,reps), main="Null distribution of KW statistic")
    abline(v=KW, lty=2 ,col="red")
    ```

!!! note "Adjustment for ties"

    When data are tied, we adjust the ranks using the **mid-ranks** for tied data as we did for the Wilcoxon test. The permutation method can also be applied to the $KW$ statistic.

    However, in order to maintain the chi-square approximation, the $KW$ statistic should be modified:

    $$
    KW_{ties} = \frac{KW}{1-\frac{\sum^g_{i=1}{(t_i^3-t_i)}}{N^3-N}}
    $$

    !!! abstract "Code"
        - If no ties, use `kruskal.test` function in the base `R` installation.
        - If ties exist, use `kruskal_test` function in the `R` add-on package `coin`.

## Multiple comparisons

We want more details of the differences to know the location of the differences.

Either we do 2-sample tests for each pair, but this will induce higher type I error. (If independent, then $1-0.95^3$)

!!! success "Aim: Control the type-I error while performing pairwise tests"

### Bonferroni Adjustment (union bound)

Use $\alpha_i = \frac{\alpha}{k(k-1)/2}$.

Bonferroni inequality (union bound): $P(E_1 \cup \cdots \cup E_m) \leq P(E_1) + \cdots + P(E_m)$

- Too conservative, high rate of type II error (false negative)
- Very wide confidence interval

### Fisher's Protected Least Significant Difference (LSD) (use MSE instead of two-sample error in t-test)

1. Perform F-test for equality of means (ANOVA). 
2. If result of F-test is to accept H0 (all same), then stop.
3. Perform all pairwise t-tests at significance level $\alpha$.

$$
|\bar{X_i} - \bar{X_j}| \geq t_{\alpha/2, N-k} \sqrt{MSE(\frac{1}{n_i} + \frac{1}{n_j})}
$$

- Can use $\alpha/2$ significance level
- But type I error may be greater than $\alpha$


!!! note "Why use t?"

    - Difference follows normal distribution
    - MSE is chi squared.
    - $\frac{Z}{\sqrt{\chi^2/n}} \tilde t_n$

#### Rank-based Fisher's LSD Procedure for Large Samples

- Benefit of using ranks: Less sensitive to outliers

1. Perform KW-test for equality
2. Stop if result is to accept H0
3. Perform all pariwise $z$-tests at significance level $\alpha$.

$$
|\bar{R_i} - \bar{R_j}| \geq z_{\alpha/2} \sqrt{S_R^2(\frac{1}{n_i} + \frac{1}{n_j})}
$$

- Can use normal critical value because 
- $S^2_R = \frac{N(N+1)}{12}$ or $\frac{N(N+1)}{12} - \frac{\sum^g_{i=1}(t_i^3-t_i)}{12(N-1)}$ (if ties exist)
- $Var(\bar{R_i} - \bar{R_j}) = S_R^2(\frac{1}{n_i} + \frac{1}{n_j})$

!!! note "Why use normal instead of t?"

    - $S_R^2$ is constant (variance of ranks is constant)
    - The only random thing is $|\bar{R_i} - \bar{R_j}|$, which is normal (large sample approximation).

### Tukey's Honesty Significant Difference (HSD) (largest difference between sample means among all pairs)

!!! warning "Suppose populations are normally distributed and sample sizes are equal."

$$
Q = \max_{i=j}\frac{\sqrt{n}|\bar{X_i} - \bar{X_j}|}{\sqrt{MSE}}
$$

MSE is for combined sample.

$$
|\bar{X_i} - \bar{X_j}| \geq q(\alpha, k, df)\sqrt{\frac{MSE}{n}}, df = N-k
$$

- Look up $q$-distribution table
- Guarntee experiment-wise error rate to be exactly $\alpha$.

??? note "Unequal sample size: Tukey-Kramer procedure"
    $$
    |\bar{X_i} - \bar{X_j}| \geq q(\alpha, k, df)\sqrt{\frac{MSE}{2}(\frac{1}{n_i} + \frac{1}{n_j})}
    $$

    Conservative, error rate $\leq \alpha$.

#### Rank-based Tukey's HSD Procedure for Large Samples

$$
|\bar{R_i} - \bar{R_j}| \geq q(\alpha, k, \infty)\sqrt{\frac{S^2_R}{n}}
$$

$df = \infty$ because of large samples. 

??? note "Unequal sample size"
    $$
    |\bar{R_i} - \bar{R_j}| \geq q(\alpha, k, df)\sqrt{\frac{S^2_R}{2}(\frac{1}{n_i} + \frac{1}{n_j})}
    $$

### Tukey's HSD vs Fisher's LSD

- HSD: $H_0$: All pairs are the same
- LSD: $H_0$: All pairs are the same, $\mu_1 = \mu_2 = \cdots$.

The $H_0$ for LSD spans over a greater space. So the Type I error of LSD will be smaller(?)

Two-step HSD: Evaluate the significance of the one-way ANOVA before proceeding to HSD. Replace $q(\alpha, k, df)$ with $q(\alpha, k-1, df)$.

## Multiple comparison permutation tests

Bonferroni Permutation Tests:

- Perform 2-sample permutation tests (e.g. permutation $t$-test or Wilcoxon rank sum test) on all pairs of treatments and compare the permutation $p$-value for each pair with the significance level $\alpha' = \frac{\alpha}{k(k-1)/2}$
    - Note that when comparing treatment $i$ to treatment $j$, should re-rank from $[1, n_i+n_j]$.
- If the $p$-value for one pair is less than $\alpha'$, declare significant difference between this pair of treatments.

Fisher's Protected LSD Permutation Tests:

- Perform a permutation test to test if there is any difference among the $k$ treatments at significance level $\alpha$.
- Obtain all (or a random sample) of permutations of the data from each pair of treatments. For each permutation, compute $T_{ij}$ for each $ij$ pair.
- $p$-value for comparing treatment $i, j$ is the fraction of permutations for which $|T_{ij}| \geq T_{ij, obs}$

Tukey's HSD Permutation test

- Obtain all (or a random samplle) of permutation of the data as in the permutation F-test. For each permutation, compute Q.
- Declare {ij} to e significantly different if $Q_{obs} \geq q*(\alpha)$

## Ordered alternatives: JT test

$H_1: F_1(x) \geq F_2(x) \geq \cdots \geq F_k(x)$

JT statistic: $T=\sum_{i \le j} T_{ij}$, where $T_{ij}$ is Mann-Whitney statistic.

Larger $T$ suppors $H1$.

JT test: Permutation test based on JT statistic:

- Compute $JT_{obs}$ based on the observed data
- Find all (or a random sample) of the permutations of the data and compute $JT$ for each permutation
- Compute the $p$-value as the fraction of the permutations for which $JT \geq JT_{obs}$

### Large sample approximation

Use the expected value and variance of the Mann-Whitney statistic...

TODO
