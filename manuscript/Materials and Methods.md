
### Data collection






### Model simulations




Allometrics are expressed as a power function: 
$$ y = aX^b $$
For parameter estimation, both predict and response variables are log-transformed to linearize the relationship. 
$$ \log (y) = \log(a) + b \log (X) + \varepsilon $$
where the error term follows a normal distribution, $\varepsilon \sim N(0,\sigma^2)$. This implies that $y$ follows a log-normal distribution on the original scale. To obtain predictions for $y$, the model is back-transformed by exponentiation. Due to the log-normal error structure, however, a bias correction is required. The expected values of $y$ is therefore 
$$ E[y] = \exp(\log(a) +\log(X) + \tfrac{1}{2}\sigma^2)$$
which can also be written as:
$$
y = aX^b \cdot \exp\left(\tfrac{1}{2}\sigma^2\right)
$$
This correction accounts for the bias introduced when transforming predictions back from the logarithmic scale.