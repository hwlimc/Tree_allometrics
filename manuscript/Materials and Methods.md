
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




### Stand Density Index

The stand density index (SDI) was calculated following Reineke (1933):

$$
SDI_i = N_i \left(\frac{D_i}{D_0}\right)^{\beta}
$$

where \(N_i\) is stand density (trees ha\(^{-1}\)),
\(D_i\) is quadratic mean diameter,
\(D_0\) is a reference diameter (25.4 cm),
and \(\beta\) is the self-thinning exponent.

### Species-specific maximum density

For each species (or functional group) \(g\), the maximum stand density index was estimated as:

$$
SDI_{\max,g}
=
Q_{0.99}
\left(
SDI \mid g
\right)
$$

where \(Q_{0.99}\) denotes the 99th percentile of observed SDI values.

### Relative Stand Density

Relative stand density was calculated as:

$$
RSD_i
=
\frac{SDI_i}
     {SDI_{\max,g}}
$$

Values of \(RSD = 1\) indicate the species-specific self-thinning boundary, whereas values below 1 indicate stands below the maximum attainable density.