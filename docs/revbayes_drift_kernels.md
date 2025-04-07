# Overview of RevBayes driftkernels
## Continuous variables

- $\lambda_{ij}$ for $i, j \in \{0, 1, 2\}$ (`switch_rates_pos`): $\mathrm{Dirichlet}(1,1,1,1)$ prior with drift kernel $\mathrm{Dirichlet}(\alpha x)$ where $x$ is the current value and $\alpha$ is scaling parameter.
  - $\alpha = 10$ with weight $2$
  - $\alpha = 25$ with weight $5$
- $\beta$ (`phy_scale`): $\exp(1)$ prior with drift kernel $\mathrm{Reciprocal}(x\exp(-\lambda / 2), x  \exp(\lambda / 2))$ where $x$ is current value and $\lambda$ is a scaling parameter
  - $\lambda$ is default value (presumably 1?) with weight $1$. Asked on RevBayes slack
- $\mu$ (`clock_host`): $\exp(10)$ prior with drift kernel $\mathrm{Reciprocal}(x\exp(-\lambda / 2), x  \exp(\lambda / 2))$ where $x$ is current value and $\lambda$ is a scaling parameter
  - $\lambda$ is default value (presumably 1?) with weight $2$. Asked on RevBayes slack
  - $\lambda = 0.2$ with weight $5$.
