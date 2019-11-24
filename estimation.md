# Parameter Estimation

## Distances Between Distributions

1. __Total Variation Distance__: The total variation distance between distributions $P$ and $Q$, with density functions $p$ and $q$, is defined to be
$$ TV(P,Q) = \sup_{A\subset E}\lvert\mathbb{P}_p(A) - \mathbb{P}_q(A)\vert. $$
This is equivalent to half the $L^1$ distance between $p$ and $q$:
$$ TV(P,Q) = \frac{1}{2}\int_E\lvert p(x) - q(x)\rvert dx. $$
    1. This is a genuine metric.
    2. Unfortunately, it is hard to estimate.
2. __Kullback-Leibler Divergence__: The Kullback-Leibler (known as relative entropy in Information Theory) divergence between $P$ and $Q$ is defined to be $$ KL(P\Vert Q) = \int_Ep(x)\log\Big(\frac{p(x)}{q(x)}\Big)dx, $$
where we assign the value $\infty$ if the support of $p$ is _not_ contained in the support of $q$ (if it is, then anywhere $q=0$, we will also have $p=0$ and thus the points at which the integrand is not defined will all be removable discontinuities).
    1. While positive semi-definite, KL-divergence is not a true metric, since it is anti-symmetric.  It also fails to satify a trangle inequality.
    2. It is, however, an expectation. Hence, it can be replaced with a sample mean and estimated.
        * Professor Rigollet calls the act of replating an expectation with a sample mean (i.e., the application of LLN) "the statistical hammer."  The implication here is that it's our simplest (and often only) tool.

### Examples

1. Let $X_n = \text{Poi}(1/n)$ and let $\delta_0$ be a point mass centered at 0.  Then $TV(X_n,\delta_0) \to 0$.
2. Let $P = \text{Bin}(n,p)$, $Q = \text{Bin}(n,q)$, where $p,q\in(0,1)$, and write their densities with one function
$$ f(p, k) = {n \choose k}p^k(1-p)^{n-k}, $$
and similarly for $f(n,q)$.  Then it is actually a pretty straightforward calculation to show that
$$ KL(P\Vert Q) = np \cdot \log\Big(\frac{p}{q}\Big)+ (n-np)\cdot\log\Big(\frac{1-p}{1-q}\Big). $$
3. Let $P = N(a,1)$ and let $Q = N(b,1)$.  Then (also pretty straightforward to calculate):
$$ KL(P\Vert Q) = \frac{1}{2}(a-b)^2. $$

## Maximum Likelihood Estimation

### Definitions

1. Let $X_1,X_2,\ldots,X_n$ be an _iid_ sample from a distribution with density $f(x; \theta)$.  The _Likelihood_ of the sample is $$ L(X_1,X_2,\ldots,X_n) = \prod_{i=1}^nf(X_i; \theta). $$
2. The _log-likelihood_ function, denoted $\ell(\theta)$ is $$ \ell(\theta) = \log(L(X_1, X_2, \ldots, X_n; \theta)). $$ Note we write $\ell$ as a random function of $\theta$.
3. The _Fisher Information_ is defined to be
$$ I(\theta) = E\big[\nabla\ell(\theta)(\nabla\ell(\theta))^T\big] - E\big[\nabla\ell(\theta)\big]E\big[\nabla\ell(\theta)\big]^T = -E\big[\mathbf{H}\ell(\theta)\big], $$
where in this case the likelihood is of a one-element sample, and the bold __H__ denotes the Hessian operator.  In one dimension, this reduces to
$$ I(\theta) = -E(\ell''(\theta)).$$
Equivalently, we also have $$ I(\theta) = Var(\ell'(\theta)). $$ This latter definition is usually harder to work with, but has a more direct connection to maximum likelihood estimators.

Throughout, we will be discussing ways to estimate the value of a "true" parameter $\theta^\ast$ of a distribution $\mathbb{P}_{\theta^\ast}$, given a model $(E, \{\mathbb{P}_\theta:\theta\in\Theta\})$.  A noble goal might be to build an estimator $\widehat{TV}(P_\theta, P_{\theta^\ast})$ and compute the argmin using this estimator.  However, $TV$ distance is hard to estimate in general, so we use $KL$-divergence instead.  Since this function is an expectation, it can be replaced by a sample mean (using LLN), and is therefore easy to estimate.

For the rest of this section, suppose we are estimating a distribution $\mathbb{P} = \mathbb{P}_{\theta^\ast}$ with a parametric family of distributions $\{\mathbb{P}_\theta : \theta\in\Theta\}$.  We will proceed to do this by estimating the minimizer (argmin) of $KL(\mathbb{P}_\theta, \mathbb{P})$, which is $\theta^\ast$, by the positive semidefiniteness (or nonnegative definiteness?) of $KL$.

The strategy for doing so will involve first estimating KL divergence and finding the minimizer of that estimator $\widehat{KL}$.  That the argmin of $\widehat{KL}$ converges to the argmin of $KL$ follows from "nice analytic properties" of these functions.  I'm guessing that $KL$ is at least $C^1$ and the convergence is relatively strong.

### Estimating $KL$ Divergence

Recall that $KL(\mathbb{P}_\theta, \mathbb{P})$ is an expectation: if $f_\theta$ and $f$ are the densities of $\mathbb{P}_\theta$ and $\mathbb{P}$, respectively, then
$$ KL(\mathbb{P}, \mathbb{P}_\theta) = E_{\mathbb{P}}\bigg(\log\bigg(\frac{f(x)}{f_\theta(x)}\bigg)\bigg) = E_{\mathbb{P}}(\log(f(x))) - E_{\mathbb{P}}(\log(f_\theta(x))). $$
As a function $\theta\mapsto KL(f_\theta, f)$, this has the form
$$ KL(\mathbb{P}, \mathbb{P}_\theta) = \text{"constant"} - E_{\mathbb{P}}(\log(f_\theta(x))). $$
Thus, by LLN, we have
$$ \widehat{KL}(\mathbb{P}, \mathbb{P}_\theta) = \text{"constant"} - \frac{1}{n}\sum_{i=1}^n\log(f_\theta(X_i)). $$

### Finding the Minimum of $\widehat{KL}$

Starting with the above equation, we have
$$\begin{aligned}
\min_{\theta\in\Theta}\widehat{KL}(\mathbb{P}, \mathbb{P}_\theta) &\Leftrightarrow \min_{\theta\in\Theta}- \frac{1}{n}\sum_{i=1}^n\log(f_\theta(X_i)) \\
&\Leftrightarrow \max_{\theta\in\Theta}\frac{1}{n}\sum_{i=1}^n\log(f_\theta(X_i)) \\
&\Leftrightarrow \max_{\theta\in\Theta}\log\bigg(\prod_{i=1}^n p_\theta(X_i)\bigg) \\
&\Leftrightarrow \max_{\theta\in\Theta}\prod_{i=1}^n p_\theta(X_i).
\end{aligned}$$
Therefore, the minimizer of $\widehat{KL}$ is the maximum likelihood estimator $\hat{\theta}$ of $\theta^\ast$.  Furthermore (avoiding a bunch of details necessary for this implication), we have
$$ \hat{\theta} \overset{(p)}{\to} \theta^\ast. $$

### Examples of Maximum Likelihood Estimators

1. __TODO__

### The Asymptotic Variance of MLE

The MLE is not only consistent, but also satisfies a central limit theorem:
$$ \sqrt{n}(\hat{\theta} - \theta^\ast) \overset{(d)}{\to} N\big(0, V(\theta^\ast)\big), $$
where $V(\theta^\ast)$ represents the asymptotic variance of $\hat{\theta}$.  But what is this asymptotic variance?!?  Turns out that under some mild conditions, the asymptotic variance of $\hat{\theta}$ is known.

__Theorem__ Assume the following.
1. $\theta^\ast$ is identifiable.
2. $\theta^\ast$ is an interior point of $\Theta$.
3. The Fisher information matrix $I(\theta)$ is invertible in a neighborhood of $\theta^\ast$.
4. All the functions involved are "nice".
5. The support of $\mathbb{P}_\theta$ does not depend on $\theta$.

Then $$ V(\hat{\theta}) = I(\theta^\ast)^{-1}. $$

__Proof__
Write $\ell_i(\theta) = \log f_\theta(X_i)$.  We start with a couple of observations:
1. Since $\hat{\theta}$ is the unique maximizer of $\log(L_n(X_1,X_2, \ldots,X_n; \theta))$, $$ \frac{d}{d\theta}\bigg\lvert_{\theta=\hat{\theta}}\sum_{i=1}^n\ell_i = \sum_{i=1}^n\ell_i(\hat{\theta}) = 0. $$
2. Since $\theta^\ast$ is the unique minimizer of $KL(\mathbb{P}_\theta, \mathbb{P})$ and this differs from $E(\ell(\theta))$ by a constant, we have $$ E(\ell'(\theta^\ast)) = \frac{d}{d\theta}\bigg\lvert_{\theta=\theta^\ast}E(\ell(\theta)) = 0. $$

Now, we start with a Taylor expansion at $\theta^\ast$:
$$ 0 = \sum_{i=1}^n\ell_i'(\hat{\theta}) = \sum_{i=1}^n\big[ \ell_i'(\theta^\ast) + (\hat{\theta} - \theta^\ast)\ell_i''(\theta^\ast) + \cdots \big]. $$
Therefore scaling and applying observation 1, we have
$$\begin{aligned}
0 &= \frac{1}{\sqrt{n}}\sum_{i=1}^n\big[ \big(\ell_i'(\theta^\ast) - E[\ell_i'(\theta^\ast)]\big) + (\hat{\theta} - \theta^\ast)\ell_i''(\theta^\ast) + \cdots \big] \\
&\approx \frac{1}{\sqrt{n}}\sum_{i=1}^n\big(\ell_i'(\theta^\ast) - E[\ell_i'(\theta^\ast)]\big) + \frac{1}{\sqrt{n}}\sum_{i=1}^n(\hat{\theta} - \theta^\ast)\ell_i''(\theta^\ast) \\
&= \frac{1}{\sqrt{n}}\sum_{i=1}^n\big(\ell_i'(\theta^\ast) - E[\ell_i'(\theta^\ast)]\big) + \sqrt{n}(\hat{\theta} - \theta^\ast)\cdot\frac{1}{n}\sum_{i=1}^n\ell_i''(\theta^\ast).
\end{aligned}$$
By CLT, the term on the left converges to $N(0, I(\theta^\ast))$.  By LLN, the term $n^{-1}\sum_i\ell''_i(\theta^\ast)$ converges to $E(\ell''(\theta^\ast)) = -I(\theta^\ast)$.  Therefore, rearranging, we have
$$ I(\theta^\ast)\cdot\sqrt{n}(\hat{\theta} - \theta^\ast) \overset{(d)}{\to} N(0,I(\theta^\ast)), $$
therefore,
$$ \sqrt{n}(\hat{\theta} - \theta^\ast) \overset{(d)}{\to} N(0,I(\theta^\ast)^{-1}). $$

__Remark:__ This proof only works in one dimension.  In higher dimensions, there is a lack of commutativity that results in a more complicated expression in the end.  (__TODO:__ Write up something about it)

__Remark:__ Recall the original definition of Fisher information as the Hessian of log-likelihood.  This adds geometric intuition to the result: If the log-likelihood is more tightly curved at $\theta^\ast$, then MLE will vary less around the maximum and vice versa.  The word "information" is also more than superficial with this in mind; i.e., more "information" iff less variance, which translates to tighter confidence intervals around MLE.

### Examples of Fisher Information

1. __TODO: Match the MLE examples above.__

## Method of Moments

* Requires model to be well-specified (unlike MLE, which will always find the nearest distribution $\mathbb{P}_\theta$ to $\mathbb{P}$).
* Computationally simpler though.
* The idea is we estimate the moments of $\mathbb{P}$ with the empirical moments $$ \widehat{m}_k = \frac{1}{n}\sum_{i=1}^nX_i^k $$
* By LLN, these converge to the moments of $\mathbb{P}$ (provided the model is well specified).

Here is how it works.  Suppose $\Theta \subset\mathbb{R}^d$ and write $$ M(\theta) = (m_1(\theta),\ldots,m_d(\theta)). $$ Assume $M$ is one-to-one, so that we can write $$ \theta = M^{-1}(m_1(\theta), \ldots, m_d(\theta)). $$  Then the _moments estimator_ is $$ \widehat{\theta}_n^{MM} = M^{-1}\big(\widehat{m}_1, \ldots, \widehat{m}_d\big). $$

We can generalize this to other functions $g_1(x), \ldots, g_d(x)$ which specify $\theta$, i.e., $$ \theta = M^{-1}(m_1(\theta), \ldots, m_d(\theta)), $$ where for each $k$, $$ m_k(\theta) = \mathbb{E}_\theta[g_k(X)]. $$  Then the _generalized method of moments estimator_ is $$ \widehat{\theta}_n^{GMM} = M^{-1}\big(\widehat{m}_1, \ldots, \widehat{m}_d\big), $$ where for each $k$, $$ \widehat{m}_k = \frac{1}{n}\sum_{i=1}^ng_k(X_i). $$

__Example:__ To see a simple example of why we might want to generalize beyond simply estimating moments directly, consider the normal distribution $N(\mu,\sigma^2)$.  The GMM estimator has $g_1(x) = x$ and $g_2(x) = x^2 - x$.

### Asymptotic Normality of GMM estimators

__Theorem:__
* Assume $M$ is one-to-one and $M^{-1}$ is continuously differentiable in a neighborhood of $\theta^\ast$.
* Let $\Sigma(\theta)$ be the covariance matrix of the vector $(g_1(X_1), \ldots, g_d(X_2)$ (assume this exists).

Then $$ \sqrt{n}\big(\widehat{\theta}_n^{GMM} - \theta^\ast\big) \overset{(d)}{\to} N\big(0, \Gamma(\theta^\ast)\big) \; \text{w.r.t. }\mathbb{P_{\theta^\ast}}, $$ where $$ \Gamma(\theta) = \bigg[\frac{\partial M^{-1}}{\partial \theta}\big(M(\theta)\big)\bigg]^T\Sigma(\theta)\bigg[\frac{\partial M^{-1}}{\partial \theta}\big(M(\theta)\big)\bigg]. $$

### MLE versus GMM

* In general, the MLE is more accurate.
* MLE still gives good results if model is misspecified
* Computational issues: Sometimes, the MLE is intractable but MM is easier (polynomial equations)

## M-Estimation

Suppose we are agnostic of any statistical model, and/or the quantity we are most interested in estimating is _not_ simply the parameter of a distribution.  In this case, we can still estimate the quantity by optimizing a suitable objective (e.g., minimizing a cost function).  This is called _M-Estimation_ (the M stands for maximum or minimum), and is the framework for "traditional" (not statistically motivated) machine learning.  The framework is as follows.
1. Let $X_1, X_2, \ldots, X_n$ be an _iid_ sample from an unspecified probability distribution $\mathbb{P}$.
2. Let $\mu^\ast$ be some parameter associated to the sample, e.g., some summary statistic.
3. Find a function $\rho:E\times \mathcal{M}  \to \mathbb{R}$, where $\mathcal{M}$ is the set of all possible values for $\mu$, such that the function $$ Q(\mu) = \mathbb{E}(\rho(X_1, \mu)) $$ achieves its minimum (or maximum) at $\mu^\ast$.
4. Replace the estimation with an average and proceed as with MLE.

Examples:
1. Let $E = \mathcal{M} = \mathbb{R}^d$ and let $\mu^\ast = \mathbb{E}(X)$.  An M-estimator is $\rho(x,\mu) = \lVert x - \mu \rVert_2^2$.
2. Let $E = \mathcal{M} = \mathbb{R}^d$ and let $\mu^\ast$ be a median of $\mathbb{P}$.  An M-estimator is $\rho(x,\mu) = \lVert x - \mu \rVert_1^1$.
3. Let $E = \mathcal{M} = \mathbb{R}$ and let $\mu^\ast$ be the $\alpha$-quantile of $\mathbb{P}$.  Then an M-estimator is $\rho(x, \mu) = C_\alpha(x-\mu)$, where $$ C_\alpha(x) = \left\{\begin{matrix}-(1-\alpha)x & : & x < 0\\ \alpha x & : & x \geq 0\end{matrix}\right., $$ This function is called a _check function_.
>

### Asymptotic Normality of M-estimators

In the case of MLE, we have asymptotic normality and known asymptotic variance (inverse fisher information).  To what extent do these properties generalize to M estimators?  It turns out they generalize quite well.  We will have asymptotic normality for M-estimators, and the asymptotic variance will have an expression only marginally less concise than that of the MLE (this is probably subject to some smoothness conditions on $\rho$).  First, we make the following definitions.  In one dimension, let $$ J(\mu) = \frac{\partial^2}{\partial\mu\partial\mu^T}Q(\mu) = \mathbb{E}\bigg[\frac{\partial^2}{\partial\mu\partial\mu^T}\rho(X_1,\mu)\bigg] $$ and let $$ K(\mu) = \text{Cov}\bigg[\frac{\partial}{\partial\mu}\rho(X_1,\mu)\bigg]. $$  In higher dimensions, $$ J(\mu) = \mathbb{E}[\mathbf{H}\rho] $$ is the expected curvature of loss and $$ K(\mu) = \text{Cov}[\nabla_\mu\rho(X_1,\mu)] $$ is the covariance matrix of the loss gradient (as a function of $\mu$ only).

__Remark:__ In the case of MLE, $J(\theta) = K(\theta) = I(\theta)$.

__Theorem:__ With notation as above, assume the following.
1. $\mu^\ast$ is the unique minimizer of $Q$;
2. $J(\mu) is invertible in a neighborhood of $\mu^\ast$;
3. A "few more technical conditions." (e.g., twice-differentiability of $\rho$, inverse of $J$ is continuous, etc.).
Then $\widehat{\mu}_n$ satisfies $$ \widehat{\mu}_n \overset{(p)}{\to} \mu^\ast $$ and $$ \sqrt{n}(\widehat{\mu}_n - \mu^\ast) \overset{(d)}{\to} N\big(0, J(\mu^\ast)^{-1} K(\mu^\ast) J(\mu^\ast)^{-1}\big). $$

The proof of this theorem is very similar to the MLE case in one dimension.  __TODO:__ Find a proof for higher dimensions, or at least a sketch of one.
