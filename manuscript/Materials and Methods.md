
### Data collection






### Model simulations




Allometrics are expressed as a power function: $y = aX^b$. For parameter estimation, both predict and response variables are log-transformed to linearize the relationship. 

$$
\log (y) = \log(a) + b \log (X) + \varepsilon
$$

, where the error term follows a normal distribution, $\varepsilon \sim N(0,\sigma^2)$. This implies that $y$ follows a log-normal distribution on the original scale. To obtain predictions for $y$, the model is back-transformed by exponentiation. Due to the log-normal error structure, however, a bias correction is required. The expected values of $y$ is therefore, $E[y] = \exp(\log(a) +\log(X) + \tfrac{1}{2}\sigma^2)$, which can also be written as:
$y = aX^b \cdot \exp\left(\tfrac{1}{2}\sigma^2\right)$. This correction accounts for the bias introduced when transforming predictions back from the logarithmic scale.

### Stand Density Index

The stand density index (SDI) was calculated following Reineke (1933):

$$
SDI_i = N_i \left(\frac{D_i}{D_0}\right)^{\beta}
$$

where $N_i$ is stand density (trees ha$^{-1})$, $D_i$ is quadratic mean diameter, $D_0$ is a reference diameter (25.4 cm), and $\beta$ is the self-thinning exponent.

### Species-specific maximum density
For each species (or functional group) (g), the maximum stand density index was estimated as the 99th percentile of observed SDI values for each species.

### Relative Stand Density
Relative stand density was calculated as:

$$
RSD_i = \frac{SDI_i}{SDI_{\max,g}}
$$

Values of \(RSD = 1\) indicate the species-specific self-thinning boundary, whereas values below 1 indicate stands below the maximum attainable density.

## Bayesian inference

Our data have a hierarchical structure and follow non-linear patterns, making robust paramerization challenging despite the relatively large sample size. To obstain robust parameter estimates and evaluate model performance, while effectively describing our data structure, we performed Bayesian inference using the brms package in R. Model parameters were estimated from posterior distributions obtained via Markov chain Monte Carlo (MCMC) sampling using four chains with 4000 iterations each.

First, MCMC convergence and sampling quality were assessed by i) checking divergent transitions to identify potential failures of the Hamiltonian Monte Carlo sampler to adequately explore the posterior distribution.by examining the potential scale reduction factor to assess whether multiple MCMC chains converged to the same posterior distribution; bulk and tail effective sample size to assess the effective number of independent samples available for estimating the central and tail regions of the posterior distributions, respectively. 

Second, posterior draws were extracted to characterize the uncertainty and variability of the estimated model parameters. Using the posterior draws, we further evaluated model fit and potential biases by comparing posterior predictive distributions with the observed data.

Finally, models with different hierarchical random-effects structures were compared using Pareto-smoothed importance sampling leave-one-out cross-validation implemented in the **loo** package. This approach estimates how well each model is expected to predict new observations by evaluating its ability to predict each observed data point after excluding that observation from model fitting. Models were then ranked according to their predictive performance, and the model with the best predictive performance was selected as the final model.
