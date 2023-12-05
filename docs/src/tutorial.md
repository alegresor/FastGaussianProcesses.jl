```@meta
DocTestSetup = quote
    using FastGaussianProcesses
    using QMCGenerators
    using Printf
    Base.show(io::IO, f::Float64) = @printf(io, "%1.3f", f)
end
```

# Tutorial

```@contents
Pages = ["tutorial.md"]
Depth = 4
```

To begin, install this package with 

```julia 
Pkg.add("FastGaussianProcesses")
```

and then import via

```julia
using FastGaussianProcesses
```

To construct specific sampling sequences we also need to be 

```julia 
using QMCGenerators
```

## Plots

```jldoctest plots; output = false
PLOTDIR = joinpath(@__DIR__,"src/assets")
f1s(x::Vector{Float64}) = [x[1]*sin(10*x[1])]
# output
f1s (generic function with 1 method)
```

### 1D Plot 

```jldoctest plots; output = false
gp = FastGaussianProcess(f1s,RandomShift(LatticeSeqB2(1),1,7),2^2;verbose=0)
plot_gp_optimization(gp,figpath=joinpath(PLOTDIR,"gp.1s.latticeseqb2.svg"),px_per_unit=4)
# output
Figure()
```

![image](./assets/gp.1s.latticeseqb2.svg)
