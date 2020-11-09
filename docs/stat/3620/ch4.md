# Chapter 4: Paired Comparisons and Block Designs

### Wilcoxon signed-rank test for paired samples

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

### Large sample approximation of Wilcoxon signed-rank statistic

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

## Page’s Test for Ordered Alternatives for RCB Designs

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

### Cochran’s Q Test (Friedman’s test with ties for 0/1 responses)

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

###  Kendall’s W Test (Test for agreements among judges)

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