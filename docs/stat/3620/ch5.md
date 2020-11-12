# Chapter 5: Tests for Trends and Association

## Pearson relation coefficient

**Pearson product moment correlation**

![Pearson product moment correlation](https://upload.wikimedia.org/wikipedia/commons/thumb/d/d4/Correlation_examples2.svg/1200px-Correlation_examples2.svg.png)

### Parametric test for zero correlation

If $(X_i, Y_i)$ is a random sample froma  bivariate normal distribution, the test statistic for the null hypothesis $H_0: \phi=0$ is 

$$
t_{corr} = \sqrt{\frac{n-2}{1-r^2}}r \tilde t(n-2)
$$

Proof comes from test for zero slope of $Y = \alpha + \beta X + \epsilon$, the test statistic for correlation is the same as the test statistic for zero slope. $t_{corr} = t_{slope}, \beta = 0$

### Permutation test for zero slope

