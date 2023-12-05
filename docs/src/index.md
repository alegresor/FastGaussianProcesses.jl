# Fast Gaussian Processes

See the [Tutorial](@ref) for instructions on how to use this package. 

`FastGaussianProcesses.jl` implements fast construction of Gaussian Process surrogates by matching shift invariant kernels to quasi-random (low discrepancy) sampling locations. For a standard GP fit to $n$ points, computing the posterior mean and variance costs $\mathcal{O}(n^3)$. For a Fast GP, these cost only $\mathcal{O}(n \log n)$. Note that Fast GPs require control of sampling locations.

This package is compatible with two flavors of quasi-random sequences: *Lattice rules* and *digital sequences*, both implemented in base 2. Sequence generators are implemented in [`QMCGenerators.jl`](https://github.com/alegresor/QMCGenerators.jl). 

## Example

Below $m_n$ is the posterior mean, $\sigma_n^2$ the posterior variance, and  

$$f(x) = ???.$$

![image](./assets/gp.1s.latticeseqb2.svg)

## References

1. Rathinavel, J. (2019). [Fast automatic Bayesian cubature using matching kernels and designs.](https://www.proquest.com/openview/0497dd427a10e6f26123b7cf5f2960f4/1?casa_token=JZfU3Big1e8AAAAA:mHuoTtMTvWnTT8q6l3AoR-eTzbbXUVAFRzSLHK8JAdQh2vvKqG3XcI_j1VV7BLKYyT43aUavfv1Q&cbl=18750&diss=y&pq-origsite=gscholar&parentSessionId=b%2BkcutbmKeHFmJujGPqXFVjhQme%2FYaRlnpP6nqlcFAk%3D) Illinois Institute of Technology.

2. Jagadeeswaran, R., & Hickernell, F. J. (2019). [Fast automatic Bayesian cubature using lattice sampling.](https://link.springer.com/article/10.1007/s11222-019-09895-9) Statistics and Computing, 29(6), 1215-1229.

3. Kaarnioja, V., Kuo, F. Y., & Sloan, I. H. (2023). [Lattice-based kernel approximation and serendipitous weights for parametric PDEs in very high dimensions.](https://arxiv.org/abs/2303.17755) arXiv preprint arXiv:2303.17755.

4. Rasmussen, C. E., & Williams, C. K. (2006). [Gaussian processes for machine learning (Vol. 1, p. 159).](https://gaussianprocess.org/gpml/) Cambridge, MA: MIT press.
