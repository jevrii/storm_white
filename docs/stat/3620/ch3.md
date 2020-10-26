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

## Kruskal-Wallis Test

!!! important "Objective"
    Obtaining a statistic based on **ranks**, similar to Wilcoxcon test.

    Approximating permutation test using chi-square, useful for large sample size.

### Kruskal-Wallis statistic

Rank the data combining all treatments. Let $\bar{R_i}$ denote the mean of the ranks of treatment $i$'s observations, assuming the data have no ties.

The **Kruskal-Wallis statistic** is defined as:

$$
KW = \frac{12}{N(N+1)}\sum^k_{i=1}n_i\left(\bar{R_i}-\frac{N+1}{2}\right)^2 = \frac{12}{N(N+1)} \sum^k_{i=1} n_i (\bar{R_i})^2 - 3(N+1)
$$

- where $\frac{N+1}{2}$ is the overall mean of all ranks and $\sum^k_{i=1}n_i\left(\bar{R_i}-\frac{N+1}{2}\right)^2$ can be viewed as the $SSB$ on ranks. (recall $SSB = \sum^k_{i=1}n_i(\bar{X_i} - \bar{X})^2$)
- $\frac{12}{N(N+1)}$ is a constant scaling factor to make the $KW$ statistic to follow **approximately** a chi-square distribution under $H_0$, i.e. $KW \sim \chi^2(k-1)$.

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

### Adjustment for ties

