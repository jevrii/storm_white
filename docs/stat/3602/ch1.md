# Chapter 1: Decision Problem: Frequentist Approach

## General setup

- Parameter space: $\Theta$
    - The true parameter is some **unknown** $\theta\in\Theta$
- Sample space $S$
    - Collection of all the possible realizations $x$ of a random vector $X$
- Statistical model $f(\cdot | \theta)$
    - "True" probability function of the random vector $X$
- Action space $A$
- Loss function $L(\theta, a)$

## Frequentist approach

- Decision rule: $d(x): S \to A$
    - Action $A$ to be taken when $X$ is observed
- Risk function: $R(\theta, d) = \mathbb{E}_\theta[L(\theta, d(X))] = \int_S{L(\theta, d(x)) f(x|\theta) dx}$
    - Characterize the performance of rule $d$ under each possible $\theta\in\Theta$

!!! success "Goal: Choose one rule $d(\cdot)$ (by their risk functions) that is the "best" (in **some** sense)"

!!! example "Example: Quality control"
    From a batch of 10 batteries, draw 1 at random and test it. Either it is defective or OK. This is the only observation upon which our decision about the whole batch has to rely.

    - Sample space: $S = {0, 1}$
    - Parameter space: $\theta\in\Theta = {0, 1, \cdots, 10}$
    - Action space: $a_0$ sell all, $a_1$ scrap all
    - Loss function: $L(\theta, a_0) = 30\theta - 150$, $L(\theta, a_1) = 10$
    - Only 4 possible decision rules: $d_{00}, d_{01}, d_{10}, d_{11}$
        - Can plot the risk function for each $d$. x axis is $\theta$, y axis is the risk.
        - $R(\theta, d) = \Pr_\theta(X = 0 | \theta)L(\theta, d_1(0)) + \Pr_\theta(X = 1 | \theta)L(\theta, d_1(1))$

## Criteria for selecting decision rules

### Admissible

!!! abstract "Admissible: Not being strictly domintated by another rule"
    A rule $d$ strictly dominates another rule $d^*$ if
    
    - $R(\theta, d) \leq R(\theta, d^*)$ for all $\theta\in\Theta$, and
    - $R(\theta', d) < R(\theta', d^*)$ for some $\theta'\in\Theta$.

A rather weak (but sensible) condition on choice of decision rule... An inadmissible rule is obviously stupid! Admissibilty provides the first criterion for discarding obviously stupid rules.

### Minimax

!!! abstract "Minimaxity: Have the smallest "worst-case" risk among all the "worst-case" risks of all rules"

Formally, rule $d$ is minimax if $\sup{R(\theta, d'): \theta\in\Theta} \geq \sup{R(\theta, d): \theta\in\Theta}$

...a conservative rule?

### Unbiased

!!! abstract "Unbiased: Incurs the smallest expected loss if the loss is measured w.r.t to the true $\theta$ generating $X$."
    That is, the loss may be measured accoding to another $\theta'$. However, the probabilty distribution of the loss (that governs the observation $X$) remains the original. 

    Unbiased rule: The loss is the smallest when the two $\theta$ are the same, opposed to when the loss $\theta$ and probabilty $\theta$ are different.

Formally, a rule $d$ is unbiased if $\mathbb{E}_\theta[L(\theta', d(X))] \geq \mathbb{E}_\theta[L(\theta, d(X))] = R(\theta, d)$ for all $\theta,\theta'\in\Theta$

!!! warning
    When determining whether a rule is unbiased or not, remember to start from the definition!

!!! example
    TODO

Application: It may be intractable to find the **best** rule in the class of all rules, but we may be able to find the **best** rule in the subclass of all **unbiased** rules.

### Bayes

!!! question "If $R(\theta, d_1) > R(\theta, d_2)$ for some $\theta$, but $R(\theta, d_1) < R(\theta, d_2)$ for other $\theta$, which is better?"

Sometimes, although we do not know the true $\theta$, somehow we may have a **subjective** opinion on how the true $\theta$ is like (e.g. based on past experience). We call this the _prior knowledge_ of unknown $\theta$.

With the prior weight function $\pi(\theta)$, we can reduce the risk function to a scalar. The weight function indicates how much "belief" we have in $\theta$ being the true parameter.

Bayes risk: $r(\pi, d) = \int_\Theta{R(\theta, d)\pi(\theta)d\theta}$

Bayes rule: Rule that has the smallest Bayes risk.

!!! success "Goal: Given $\pi(\cdot)$, seek the Bayes rule (w.r.t. prior $\pi$)."

To find Bayes rules in practice, we can find $d$ that minimizes $\int_S\int_\theta L(\theta, d(x)) \pi(\theta) f(x|\theta) d\theta dx$ instead. (seems that you can just minimize the inner integral?)

??? Derivation
    TODO

## Randomized decision rule

A random mixture of rules - Choose a decision rule by random! :game_die:

> Mainly for analytical purposes... seldom put into practice

- Risk function: $R(\theta, d^*) = \sum^s_{i=1} p_i R(\theta, d_i)$
    - **Convex** combination of the risk functions of the individual rules $d_i$'s.
    - So it won't exceed the original deterministic curves
- Bayes risk: $r(\pi, d^*) = \int R(\theta, d^*) \pi(\theta) d\theta = \sum^s_{i=1}p_i\int R(\theta, d_i) \pi(\theta) d\theta = \sum^s_{i=1}p_i r(\pi, d_i) \geq \min (r(\pi, d_1), \cdots, r(\pi, d_s))$
    - Thus, the risk of a Bayes rule chosen among $(d_1, \cdots d_s)$ cannot be reduced further by making them into a randomized rule.
- Minimax: Randomized rule **can** be more minimax than deterministic minimax rules!

## Risk set

(restrict to special setting: $\Theta = {\theta_1, \theta_2}$)

2D plot: x axis $R(\theta_1, d)$, y axis $R(\theta_2, d)$.

Then each decision rule is a point in the plot. Can find the good rules easily!

All points in the convex hull -> includes both deterministic and randomized decision rules

Risk set must be convex

- Admissible rules: lower-left boundary of risk set
- Minimax: Moving the "right angle" $\max(x, y) = c$ until it touches the risk set. In some cases, draw a line $R(\theta_1, d) = R(\theta_2, d)$ and find the intersection with the lower boundary, but **this doesn't always work**.
- Unbiased rule: Nil (not using the diagram _alone_). Because the risk set only tells us the risk function, but the definition of unbiased rules requires formulae outside of the risk function.
- Bayes rule (w.r.t prior $\pi(\theta_1) = \pi_1, \pi(\theta_2) = \pi_2$): Push a line $\pi_1 R(\theta_1, d) + \pi_2 R(\theta_2, d) = c$. Minimize $c$ such that it is still touching the risk set.

!!! example "Example: Hypothesis testing"
    TODO

    Admissible rule: Lower boundary. Refer to Neyman-Pearson lemma.

    Minimax rule: $R(0, d_k) = R(1, d_k) \implies 1-\Phi(k) = \Phi(k-1) \implies \Phi(-k) = \Phi(k-1) \implies k = 0.5$

    Unbiased rule: $R(0, d_k) \leq \frac{1}{2}, R(1, d_k) \leq \frac{1}{2}$

    Randomized rule: Find the point via section formula

    Bayes rule: Minimize $\psi (1-\Phi(k)) + (1-\psi) \Phi(k-1)$. Differentiating w.r.t. $k$ gives $-\psi \phi(k) + (1-\psi) \psi(k-1) = 0$, then solve for $k$.

    - We can find the Bayes prior for a given rule. In order words, what prior needs to be exist for this rule to be a Bayes rule?
    - E.g. $k = \Phi^{-1}(0.95)$, solve for $\psi$. Then $\psi = 0.75857$
    - Size 5% test -> 75% belief on the null hypothesis. We are not making a neutral decision!
    - Equivalent to finding slope $= -\frac{\psi}{1-\psi}$ of risk set at $R(0, d) = 0.05$
