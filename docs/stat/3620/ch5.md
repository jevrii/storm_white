# Chapter 5: Tests for Trends and Association

## Pearson relation coefficient

**Pearson product moment correlation**

![Pearson product moment correlation](https://upload.wikimedia.org/wikipedia/commons/thumb/d/d4/Correlation_examples2.svg/1200px-Correlation_examples2.svg.png)

$= \pm 1$: **Perfect** linear relationship $Y = a+bX$

### Parametric test for zero correlation

If $(X_i, Y_i)$ is a random sample froma  bivariate normal distribution, the test statistic for the null hypothesis $H_0: \rho=0$ is 

$$
t_{corr} = \sqrt{\frac{n-2}{1-r^2}}r \sim t(n-2)
$$

Proof comes from test for zero slope of $Y = \alpha + \beta X + \epsilon$, the test statistic for correlation is the same as the test statistic for zero slope. $t_{corr} = t_{slope}, \beta = 0$ ($\rho$ is monotonic function of $\beta$)

(but the problem for testing zero slope is slightly different from the one for linear model, as linear model assume $X$ is fixed, while here both $X$ and $Y$ are random)

### Permutation test for zero slope

1. Calculate $\hat{\beta_{obs}}$ using the observed data
2. Keep the order of $X$ unchanged, permute the order of $Y$ observation. There are $n!$ observations.
3. Compute $\hat{\beta}$ for each possible permutation
4. Find p-value

**Large sample approximation for $r$**

$$
Z = \frac{r}{\sqrt{\frac{1}{n-1}}} = \sqrt{n-1}r \sim N(0, 1)
$$

## Spearman Rank Correlation

Issue: Pearson's correlation coefficient can only capture the linear relationship.

Solution: Consider the correlation of the ranking data instead.

**Spearman's rank correlation** is obtained by applying the Pearson correlation to the pairs $(R(X_i), R(Y_i))$.

When there is no ties in ranking, $r_s$ can be simplified in terms of the differences $D_i = R(X_i) - R(Y_i)$: $r_s = 1 - \frac{6\sum^n_{i=1}D_i^2}{n(n^2-1)}$

**Large sample approximation for $r_s$**

$$
Z = \sqrt{n-1}r_s \sim N(0, 1)
$$

## Kendall's Tau

- Concordant: $(X_i - X_j)(Y_i - Y_j) > 0$
- Discordant: $(X_i - X_j)(Y_i - Y_j) < 0$

**Kendall's tau**: $\tau = p_c = p_d$

$$
r_\tau = \frac{\sum\sum_{i<j} sgn[(X_i - X_j)(Y_i - Y_j)]}{{n\choose 2}}
$$

**Large sample approximation**

TODO

## Contingency table

Chi-square test