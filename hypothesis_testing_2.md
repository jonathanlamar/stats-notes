# Hypothesis Testing II

The other hypothesis testing doc covered some extreme basics, and left a lot of
questions unanswered (*but what about the __geometry__ of error!?*).  While I
probably won't ever have time to explore some of those concerns, this doc will
cover some more intermediate topics in hypothesis testing.

## $t$-Tests

### A Motivating Example

Suppose we are performing clinical trials, and want to measure the effectiveness
of a drug over placebo at reducing cholesterol.  Let $X_1, \ldots, X_n$ be the
test group measured drop in cholesterol and let $Y_1,\ldots, Y_m$ be the control
group measured drop in cholesterol.  Then $X_1,\ldots,X_n$ are _iid_
$N(\Delta_d, \sigma_d^2)$ and $Y_1,\ldots,Y_m$ are _iid_ $N(\Delta_c,
\sigma_c^2)$, and our null hypothesis is that $\Delta_d \geq \Delta_c$.  Using
Slutskey's Lemma and CLT, we have $$ \frac{\bar{X_n} - \bar{Y_m} - (\Delta_d -
\Delta_c)}{\sqrt{\frac{\hat{\sigma}^2_d}{n} + \frac{\hat{\sigma}^2_c}{m}}}
\overset{(d)}{\to} N(0,1). $$ _However_, this requires both sample sizes to go
to infinity, and convergence will be extremely slow if they do not do so at the
same rate.  Since we rarely have a linear relationship between control and test
group sizes (for ethical reasons), asymptotic hypothesis tests of the previous
document will not be very useful.

### Small Sample Sizes

Essentially, the above limit was a statement of the form...

### The $\chi^2$ and $t$ Distributions

Set the following notation.  Let $$ S_n = \frac{1}{n}\sum_{i=1}^n(X_i -
\bar{X}_n)^2 $$ be the sample variance of a sample of $n$ independent Gaussian
random variables with variance $\sigma^2$, and let $$ \tilde{S}_n =
\frac{1}{n-1}\sum_{i=1}^n(X_i - \bar{X}_n)^2 = \frac{n}{n-1}S_n $$ be the
unbiased estimator of sample variance.

#### The $\chi^2$ Distributions

The _$\chi^2$ distribution with $d$ degrees of freedom_, written $\chi^2_d$, is
the distribution of the sum of the squares of $d$ independent samples from
$N(0,1)$.  Equivalently (and more geometrically), if $Z \sim N(0, I_d)$ is a
$d$-dimensional Gaussian random variable with unit variance, then $\lVert
Z\rVert_2^2\sim \chi^2_d$.  This provides a geometric reason why the kurtosis
increases with $d$: The average magnitude of a Gaussian random vector will
increase with the number of dimensions (more positive numbers to sum).

Let $V \sim \chi^2_d$.  Then

* $\mathbb{E}[V] = d$; $Var(V) = 2d$.  For large $d$, we have by CLT that
* $\chi^2_d\approx N(d, 2d)$.

#### Cochrane's Theorem

__Theorem:__

1. $\bar{X}_n$ is independent of $S_n$.
2. $\frac{nS_n}{\sigma^2} \sim\chi^2_{n-1}$

#### The Student's $t$-Distribution

Let $Z$ be standard normal, let $V$ be $\chi^2_d$, and assume $Z$ and $V$ are
independent.  Then the random variable $$ t = \frac{Z}{\sqrt{V/d}} $$ has as its
distribution the _Student's $t$ distribution with $d$ degrees of freedom_.

### The Student's $t$ Test

This is the first _nonasymptotic_ test we see in this class.  One important
thing to note is that when we do nonasymptotic hypothesis testing, we cannot
escape the fact that we don't know our distribution.  This means we always place
an assumption on the underlying distribution of the sample.  In the case of the
Student's $t$-test (note the "Student" in "Student's $t$ test" is there not only
for historical reasons, but also to distinguish from Welch's $t$-test), we
assume that our data are Gaussian.

#### The One Sample Student's $t$-Test

1. Assume $X_1,\ldots,X_n$ are _iid_ Gaussian $N(\mu, \sigma^2)$.
2. The null hypothesis is that $\mu=0$, and the alternate hypothesis is either
that $\mu\neq 0$ (for the two-sided test), or $\mu\geq 0$ (for the one-sided
test).
3. The test statistic is $$ T_n := \frac{\bar{X}_n}{\sqrt{\tilde{S}_n/n}}. $$

Note that under the null, we have $$ T_n = \frac{\sqrt{n}\frac{\bar{X}_n -
\mu_0}{\sigma}}{\sqrt{\tilde{S}_n/\sigma^2}} \sim \frac{Z}{\sqrt{V}}, $$ where
$Z\sim N(0,1)$ and $V\sim\chi^2_{n-1}/(n-1)$ (by Cochrane's Theorem).  Thus, the
test statistic follows a $t$ distribution with $n-1$ degrees of freedom and
therefore its quantiles are known.

#### The Two Sample Welch's $t$-Test

Returning to out cholesterol example, we can consider the scenario of testing
for the difference of means of two samples using a $t$ distribution.  As our
null hypothesis is that $\Delta_d \geq \Delta_c$, or equivalently that
$\Delta_d - \Delta_c \geq 0$, we have a test statistic of the form $$ T_{m,n} =
    \frac{\bar{X}_n - \bar{Y}_n}{\sqrt{\frac{\hat{\sigma}^2_d}{n} +
    \frac{\hat{\sigma}^2_c}{m}}}. $$ Thus, this is a one-sided example of the
    _Welch $t$-test_.  This is opposed to the Student's $t$-test, where the test
    statistic follows a $t$ distribution.  In this case, the test statistic is
    approximately a $t$-distribution, particularly because the denominator
    involves something that is _approximately_ (and very nearly so) a $\chi^2$
    distribution.

__Theorem (Welch-Satterthwaite):__ We have $T_{m,n} \approx t_N$, where $$ N =
\frac{(\hat{\sigma}^2_d/n +
\hat{\sigma}^2_c/m)^2}{\frac{\hat{\sigma}^4_d}{n^2(n-1)} +
\frac{\hat{\sigma}^4_c}{m^2(m-1)}}\geq \min(n,m). $$

__Remark:__ If the variances are known to be equal, the test statistic becomes
exactly a $t$ distribution, hence the test becomes a two sample Student's
$t$-test.

## Tests based on MLEs

Briefly, these are some other tests.

### Wald's Test

Consider an _iid_ sample $X_1, \ldots, X_n$ with statistical model
$(E,\{\mathbb{P}_\theta\}_{\theta\in\Theta})$, where
$\Theta\subseteq\mathbb{R}^d$ and let $\theta_0$ be fixed and given.  Let
$\theta^\ast$ be the true parameter under the model.  Consider the null
hypothesis $H_0: \theta^\ast = \theta_0$ and let $\hat{\theta}_n$ be the MLE.

If $H_0$ is true, then by CLT, we have
$$
\sqrt{n}I(\theta_0)^{\frac{1}{2}}\cdot\big(\hat{\theta}_n - \theta_0)
\overset{(d)}{\to} N(0, I_d).
$$
Hence, by plugging in the MLE into the Fisher information, we have a test
statistic $T_n$ such that
$$
\underbrace{n\big(\hat{\theta}_n-\theta_0)^TI(\hat{\theta}_n)(\hat{\theta}_n-
\theta_0)}_{T_n} \overset{(d)}{\to} \chi_d^2.
$$

__Definition:__ _Wald's Test_ is any test (one or two sided) based on the above
test statistic.

### Wald's Test For Implicit Hypotheses

Similar to above, suppose our null hypothesis is of the form $H_0: g(\theta)=0$
for some continuously differentiable function $g:\mathbb{R}^d\to\mathbb{R}^k$
(with $k<d$).  Suppose an asymptotically normal estimator $\hat{\theta}_n$ is
available with asymptotic covariance $\Sigma(\theta) \in \mathbb{R}^{d\times d}$.
Let
$$
\Gamma(\theta) = \nabla g(\theta)^T\Sigma(\theta)\nabla g(\theta)
\in\mathbb{R}^{k\times k}.
$$
Then by the Delta method, we have
$$
\sqrt{n}\cdot \Gamma(\theta)^{\frac{1}{2}}\Big(g(\hat{\theta}_n) - g(\theta)\Big)
\overset{(d)}{\to} N(0, I_k).
$$
By Slutskey's Theorem, we can plug $\hat{\theta}_n$ into $\Gamma$, hence we have
a test statistic $T_n$ of the form
$$
\underbrace{n\cdot g(\hat{\theta}_n)^T\Gamma^{-1}(\hat{\theta}_n)\cdot
g(\hat{\theta})}_{T_n} \overset{(d)}{\to} \chi^2_k.
$$

__Definition:__ _Wald's Test for Implicit Hypotheses_ is any test (one or two
sided) based on the above test statistic for some function $g$.

### Likelihood Ratio Test

TODO

## Goodness of Fit Tests

### $\chi^2$ Test

TODO

TODO: Others here I barely remember from lecture.
