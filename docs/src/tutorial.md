```@meta
DocTestSetup = quote
    using FastGaussianProcesses
    using QMCGenerators
    using Printf
    using CairoMakie
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


## Kernel 1s Lines

### Lattice Kernels

```jldoctest plots; output = false
PLOTDIR = joinpath(@__DIR__,"src/assets")
plot_gp_kernel_latticeseqb2_1s_lines(figpath=joinpath(PLOTDIR,"kernels.lines.latticeseqb2.svg"))
# output
Figure()
```

![image](./assets/kernels.lines.latticeseqb2.svg)

### Digital Kernels

```jldoctest plots; output = false
plot_gp_kernel_digitalseqb2g_1s_lines(figpath=joinpath(PLOTDIR,"kernels.lines.digitalseqb2g.svg"))
# output
Figure()
```

![image](./assets/kernels.lines.digitalseqb2g.svg)


### RBF Kernels 

```jldoctest plots; output = false
plot_gp_kernel_rbf_1s_lines(figpath=joinpath(PLOTDIR,"kernels.lines.rbf.svg"))
# output
Figure()
```

![image](./assets/kernels.lines.rbf.svg)

## Kernel 1s Contours and Surfaces

### Lattice Kernels 

```jldoctest plots; output = false
PLOTDIR = joinpath(@__DIR__,"src/assets")
plot_gp_kernel_latticeseqb2_1s_contsurfs(figpath=joinpath(PLOTDIR,"kernels.contsurf.latticeseqb2.svg"))
# output
Figure()
```

![image](./assets/kernels.contsurf.latticeseqb2.svg)

### Digital Kernels

```jldoctest plots; output = false
plot_gp_kernel_digitalseqb2g_1s_contsurfs(figpath=joinpath(PLOTDIR,"kernels.contsurf.digitalseqb2g.svg"))
# output
Figure()
```

![image](./assets/kernels.contsurf.digitalseqb2g.svg)

### RBF Kernels 

```jldoctest plots; output = false
plot_gp_kernel_rbf_1s_contsurfs(figpath=joinpath(PLOTDIR,"kernels.contsurf.rbf.svg"))
# output
Figure()
```

![image](./assets/kernels.contsurf.rbf.svg)

## GPs 1s

```jldoctest plots; output = false
function f1s(x::Vector{Float64})
    [   x[1]*sin(10*x[1]), # f^{(0)}
        sin(10*x[1]) + 10*x[1]*cos(10*x[1]) # f^{(1)}
    ]
end
# output
f1s (generic function with 1 method)
```

### GP Lattice Kernel

```jldoctest plots; output = false
gp_1s_latticeseqb2 = FastGaussianProcess(f1s,RandomShift(LatticeSeqB2(1),1,7),2^2;β=[i for i=0:1,j=1:1],verbose=0)
plot_gp_optimization(gp_1s_latticeseqb2,figpath=joinpath(PLOTDIR,"optim.1s.latticeseqb2.svg"))
plot_gp_1s(gp_1s_latticeseqb2;f=f1s,β=[0,1],figpath=joinpath(PLOTDIR,"gp.1s.latticeseqb2.svg"))
# output
Figure()
```

![image](./assets/optim.1s.latticeseqb2.svg)

![image](./assets/gp.1s.latticeseqb2.svg)

### GP Digital Kernel 

```jldoctest plots; output = false
gp_1s_digitalseqb2g = FastGaussianProcess(f1s,RandomDigitalShift(DigitalSeqB2G(1),1,7),2^2;β=[i for i=0:1,j=1:1],verbose=0)
plot_gp_optimization(gp_1s_digitalseqb2g,figpath=joinpath(PLOTDIR,"optim.1s.digitalseqb2g.svg"))
plot_gp_1s(gp_1s_digitalseqb2g;f=f1s,β=[0,1],figpath=joinpath(PLOTDIR,"gp.1s.digitalseqb2g.svg"))
# output
Figure()
```

![image](./assets/optim.1s.digitalseqb2g.svg)

![image](./assets/gp.1s.digitalseqb2g.svg)

### GP RBF Kernel 

```jldoctest plots; output = false
gp_1s_rbf = GaussianProcessRBF(f1s,IIDU01Seq(1,7),2^2;β=[i for i=0:1,j=1:1],verbose=0)
plot_gp_optimization(gp_1s_rbf,figpath=joinpath(PLOTDIR,"optim.1s.rbf.svg"))
plot_gp_1s(gp_1s_rbf;f=f1s,β=[0,1],figpath=joinpath(PLOTDIR,"gp.1s.rbf.svg"))
# output
Figure()
```

![image](./assets/optim.1s.rbf.svg)

![image](./assets/gp.1s.rbf.svg)

## GPs 2s

```jldoctest plots; output = false
function f2s(x::Vector{Float64})
    [   x[1]*sin(10*x[1]) - x[2]*cos(10*x[2]), # f^{(0,0)},
        sin(10*x[1]) + 10*x[1]*cos(10*x[1]), # f^{(1,0)}
        10*x[2]*sin(10*x[2]) - cos(10*x[2]) # f^{(0,1)}
    ]
end
# output
f2s (generic function with 1 method)
```

### GP Lattice Kernel

```jldoctest plots; output = false
gp_2s_latticeseqb2 = FastGaussianProcess(f2s,RandomShift(LatticeSeqB2(2),1,7),2^8;β=[0 0; 1 0; 0 1],optim_steps=40,verbose=0)
plot_gp_optimization(gp_2s_latticeseqb2,figpath=joinpath(PLOTDIR,"optim.2s.latticeseqb2.svg"))
plot_gp_2s(gp_2s_latticeseqb2;f=f2s,β=[0 0; 1 0; 0 1],figpath=joinpath(PLOTDIR,"gp.2s.latticeseqb2.svg"))
# output
Figure()
```

![image](./assets/optim.2s.latticeseqb2.svg)

![image](./assets/gp.2s.latticeseqb2.svg)

### GP Digital Kernel

```jldoctest plots; output = false
gp_2s_digitalseqb2g = FastGaussianProcess(f2s,RandomDigitalShift(DigitalSeqB2G(2),1,7),2^8;β=[0 0; 1 0; 0 1],optim_steps=40,verbose=0)
plot_gp_optimization(gp_2s_digitalseqb2g,figpath=joinpath(PLOTDIR,"optim.2s.digitalseqb2g.svg"))
plot_gp_2s(gp_2s_digitalseqb2g;f=f2s,β=[0 0; 1 0; 0 1],figpath=joinpath(PLOTDIR,"gp.2s.digitalseqb2g.svg"))
# output
Figure()
```

![image](./assets/optim.2s.digitalseqb2g.svg)

![image](./assets/gp.2s.digitalseqb2g.svg)

### GP RBF Kernel

```jldoctest plots; output = false
gp_2s_rbf = GaussianProcessRBF(f2s,IIDU01Seq(2,7),2^8;β=[0 0; 1 0; 0 1],optim_steps=40,verbose=0)
plot_gp_optimization(gp_2s_rbf,figpath=joinpath(PLOTDIR,"optim.2s.rbf.svg"))
plot_gp_2s(gp_2s_rbf;f=f2s,β=[0 0; 1 0; 0 1],figpath=joinpath(PLOTDIR,"gp.2s.rbf.svg"))
# output
Figure()
```

![image](./assets/optim.2s.rbf.svg)

![image](./assets/gp.2s.rbf.svg)

## Logo 

In the following we must be 
```julia 
using CairoMakie
```

```jldoctest plots; output = false
m = 2; n = 2^m
x = FirstLinear(RandomDigitalShift(DigitalSeqB2G(1),1,11),m)[:,1]
kmat = Matrix{Float64}(undef,2*n,2*n)
for i1=1:n,i2=1:n,β1=0:1,β2=0:1 kmat[β1*n+i1,β2*n+i2] = ((β1+β2)==0)+kernel_digshiftinvar_s1(x[i1],x[i2],β1,β2,4) end 
fig = Figure(resolution=(400,400),backgroundcolor=:transparent)
ax = Axis(fig[1,1],backgroundcolor=:transparent)
heatmap!(ax,1:2*n,2*n:-1:1,kmat,colormap=:julia_colorscheme)#:nipy_spectral) # https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html
hlines!(ax,[n+1/2],color=:black,linewidth=2.); vlines!(ax,[n+1/2],color=:black,linewidth=2.)
hidespines!(ax); hidedecorations!(ax); hidexdecorations!(ax,grid = false); hideydecorations!(ax, ticks = false)
save(joinpath(PLOTDIR,"logo.svg"),fig)
# output
CairoMakie.Screen{IMAGE}
```

![image](./assets/logo.svg)
