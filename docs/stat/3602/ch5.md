# Chapter 5: Estimation

## Rao-Blackwell Theorem

Define $\rho^* = \rho^*(T) = \mathbb{E}[\rho(X) | T]$

1. $\rho^*$ is an estimator of $\psi(\theta)$ (it is a function of data, not depending on $\theta$)
2. $\mathbb{E}[L(\theta, \rho^*)] \leq \mathbb{E}[L(\theta, \rho)]$ (MSE)
3. If $\rho$ unbiased and $T$ complete for $\theta$, then

- $\rho^* = \rho^*(T)$ is the **only** function of $T$ which is unbiased for $\psi(\theta)$
- $\forall$ **unbiased** estimator $S$ of $\psi(\theta)$, $\mathbb{E}[L(\theta, \rho^*)] \leq \mathbb{E}[L(\theta, S)] \forall \theta$ (UMVU estimator under squared loss)

## Fisher information

Compare experiments and their likelihood functions