# Chapter 4: Paired Comparisons and Block Designs

## Paired comparisons

Paired data can be found in the following situations:

- Two measurements obtained from the same subject under two different conditions
- Two measurements obtained from matched pairs

### Parametric method: Paired t-test

Perform $t$-test on the difference data $X_i - Y_i$. $t=\frac{\bar{D}}{S_D/\sqrt{n}}$

### Paired comparison permutation test

1. Compute the sample mean of differences $D_i = X_i - Y_i$ for each pair of observed data, denoted by $\bar{D_{obs}}$.
2. Obtain the $2^n$ possible assignments (or a random sample of $R$ assignments) of plus and minus signs to the $|D_i|$'s.
3. Compute $\bar{D}$ for each possible assignment
4. Calculate the $p$-value as the fraction of $\bar{D} \geq \bar{D_{obs}}$

### Binomial (sign) test for paired samples

Apply binomial tet to the difference data $D_i$, test if population median of $D=0$. (refer to Chapter 1)

Issue: Discard a lot of information of the data, only take into acount the direction of the difference, but not the magnitude.

### Wilcoxon signed-rank test for paired samples

- Rank the absolute values of the differences $|D_i|$
- Re-attach the rank of $|D_i|$ by the sign ($+, 0, -$) of $D_i$
- Apply the permutation test on the sum of positive signed ranks, $SR_+$ or the average of signed ranks $\bar{SR}$

Assumption: $D_i$ are iid from a continuous cdf $F$ with median $\theta$, and its pdf $f$ is symmetric about $\theta$.

Signed-rank test more powerful in asymmetric case.

Can also use this tet to test for the median of a symmetric distribution based on a single sample.

??? example

    ```r
    x <- c(20, 18, 24, 14, 5, 26, 15, 29, 15, 9, 25, 31, 35, 12)
    y <- c(40, 25, 38, 27, 31, 21, 32, 38, 25, 18, 32, 28, 33, 29)

    d = x - y

    # Compute the signed ranks
    SR <- {rank(abs(d))* sign(d)}
    SRpos <- sum(SR[SR>0])

    barSR.obs = mean(SR) # average signed-rank statistic

    # Permutation test based on average signed-rank
    n <- length(SR)
        junk <- matrix(NA,2^n,n)
        for (i in 1:(2^n)) 
        {   
            remain <- i-1
            for (j in 1:n)
        {
            junk[i,j] <- remain %/% (2^(n-j)) #%/%: integer division
            remain <- remain %% (2^(n-j)) #%%: remainder
        }
    }
    signs <- junk - !junk  ### permuted signs (+/-)
    perm <- apply( t(abs(SR) * t(signs)) , 1, mean) 

    pval.twosided = mean(abs(perm) >= abs(barSR.obs))
    pval.upper = mean(perm >= barSR.obs)
    pval.lower = mean(perm <= barSR.obs)
    # pval.twosided = 2*min(pval.upper, pval.lower)
    ```

### Large sample approximation of Wilcoxon signed-rank statistic SR+

Under H0, indicator of sign of $D$ follows $Ber(0.5)$.

Then $E(SR_+) = \frac{1}{2} \sum^n_{i=1} R_i$, $Var(SR_+) = \frac{1}{4} \sum^n_{i=1} R_i^2$, $\frac{SR_+ - E(SR_+)}{\sqrt{Var(SR_+)}} \tilde N(0, 1)$

In case of ties, use average rank of tied observations.

??? example

    ```r
    x <- c(20, 18, 24, 14, 5, 26, 15, 29, 15, 9, 25, 31, 35, 12)
    y <- c(40, 25, 38, 27, 31, 21, 32, 38, 25, 18, 32, 28, 33, 29)

    d = x - y

    # Wilcoxon signed-rank test using sum of positive signed ranks
    wilcox.test(d, correct=FALSE)
    wilcox.test(x, y, paired=T, correct=FALSE)
    wilcox.test(x, y, paired=T, correct=TRUE)
    ```

## Randomized Complete Block Design (RBCD)

Apply different treatment to each person in each block. Each block contains similar people.

$X_{ij} = \mu + t_i + b_j + \epsilon_{ij}$

Test between $H_0: t_1 = \cdots = t_k$ vs $H_1$: Not all $t_i$ are the same.

When $\epsilon \tilde N(\mu, \sigma^2)$, we can carry out the $F$-test using the $F$ statistic: $F = \frac{MSB}{MSE} = \frac{SSB/(k-1)}{SSE/[(k-1)(b-1)]}$

### Permutation F-test for RCB designs

- Compute the F statistic $F_{obs}$ using the observed data
- Permute the observations within each block. There are $(k!)^b$ permutations.
- Compute the F statistic for each possible permutation.
- Calculate the p-value as the fraction of $F$'s $geq$ $F_{obs}$

Alternatively, we can also use SSB or SSX for the test statistic.

### Friedman's test for RCB design (SSB on ranks)

Test on SSB on ranks. rank within each block.

$$
FM = \frac{12b}{k(k+1)} \sum^k_{i=1}(\bar{R_{i+} - \frac{k+1}{2}})^2
$$

- For large samples, $FM \to \chi^2(k-1)$ under $H_0$.
- For $k=2$, Friedman test is equivalent to two-sided Wilcoxon signed rank test. (pair is also a block)
- Can obtain exact p-value from permutation method.

!!! note "Adjustment for ties"
    $$
    FM_{ties} = \frac{b}{\frac{1}{b}\sum^b_{j=1}S^2_{Bj}} \sum^k_{i=1}\left(\bar{R_{i+}} - \frac{k+1}{2}\right)^2
    $$

**Multiple comparisons based on Friedman Statistic for Large Samples (Re: HSD)**

$$
|\bar{R_{i+}} - \bar{R_{j+}}| \geq q(\alpha, k, \infty)\sqrt{\frac{\frac{1}{b}\sum^b_{j=1}S^2_{Bj}}{b}}
$$

## Page’s Test for Ordered Alternatives for RCB Designs

Let $R_{i+} = b\bar{R_{i+}}$ be the sum ranks for the $i$'th treatment.

$$
PG = \sum^k_{i=1} iR_{i+}
$$

Amplify the large rank sums.

- For $H_1: t_1 \leq \cdots t_k$, $PG$ tends to be large. For $H_1: t_1 \geq \cdots t_k$, $PG$ tends to be small.
- Large sample approximation: $E(PG) = \frac{bk(k+1)^2}{4}, Var(PG) = \frac{k(k^2-1)}{12}\sum^b_{j=1}S^2_{Bj}$
- Page's test will have greater power to detect difference than Friedman's test if the ordering among treatments is correct.

??? example

    ```r
    x = c(120, 208, 199, 194, 177, 195,
        207, 188, 181, 164, 155, 175,
        122, 137, 177, 177, 160, 138,
        128, 128, 160, 142, 157, 179)
    blocks = rep(1:6,4)
    grps = rep(1:4, each=6)
    k <- length(table(grps))

    # obtain the within-block ranks
    mat <- cbind(x,blocks,1:length(x))
    mat <- mat[order(mat[,2]),]
    b <- length(table(blocks))
    n <- length(x)
    for (j in 1:b) mat[n/b*(j-1)+(1:(n/b)),1] <- rank( mat[n/b*(j-1)+(1:(n/b)),1] )
    mat <- mat[order(mat[,3]),]
    Rij=mat[,1]

    # compute within-block rank sums
    junk <- table(grps)
    ranksumvec <- rep(NA,length(junk))
    for (i in 1:length(junk)) ranksumvec[i] <- sum(Rij[grps==names(junk)[i]])
    ranksumvec

    # compute observed Page's statistic
    Pageobs <- sum((1:k)*ranksumvec)

    # Permutation method based on Page's statistic
    set.seed(1234567)
    R=10000
    mat <- cbind(x,grps,blocks)
    mat <- mat[order(mat[,3]),]  # sort matrix rows by block
    k <- length(table(grps))
    b <- length(table(blocks))
    permpage <- rep(NA,R)
    for (i in 1:R)
        {
        junk <- rep(NA,b*k)
        for (j in 1:b) junk[k*(j-1)+(1:k)] <- sample(1:k,k)
        mat1 <- cbind(mat[,1],blocks,1:length(mat[,1]))
        mat1 <- mat1[order(mat1[,2]),]
        b <- length(table(blocks))
        n <- length(x)
        for (j in 1:b) mat1[n/b*(j-1)+(1:(n/b)),1] <- rank( mat1[n/b*(j-1)+(1:(n/b)),1] )
        mat1 <- mat1[order(mat1[,3]),]
        Rij=mat1[,1]	
        
        temp <- table(junk)
        sumvec <- rep(NA,length(temp))
        for (j in 1:length(temp)) sumvec[j] <- sum(Rij[junk==names(temp)[j]])
        
        permpage[i] <- sum((1:k)*sumvec)
        }

    # tests if means are arranged t1 > t2 > t3 > t4, i.e. response decreases as changing trt from 1 to 4
    pvalue = mean(permpage <= Pageobs)
    pvalue

    # correlation of treatment i with its rank sum R_i+
    cor(1:4, ranksumvec)
    ```

## Special cases of Friedman's test

### Cochran's Q Test (Friedman's test with ties for 0/1 responses)

Equivalent to Friedman's test.

??? example
    ```r
    x = c(1,0,0,1,1,1,1,0,1,1,
          1,1,1,1,0,1,1,0,0,0,
          1,1,0,1,1,1,1,1,0,1,
          1,1,0,0,0,1,1,1,0,1)
    blocks = rep(1:10,4)
    grps = rep(1:4,rep(10,4))

    # Friedman rank sum test based on chi-square approximation
    friedman.test(x, grps, blocks)  
    ```

###  Kendall's W Test (Test for agreements among judges)

Agreement: Judges are the **blocks**, interns are the **treatments**. 

In other words, the interns can be distinguished well.

??? example
    ```r
    x = c(5.5, 7.0, 8.5,
      5.0, 6.0, 6.0,
      6.0, 8.0, 9.5,
      4.5, 7.5, 8.0,
      9.0, 9.5, 9.0)
    grps = rep(1:5, each=3) 
    blocks = rep(1:3, 5) 

    # Kendall's W test
    friedman.test(x, grps, blocks)
    ```