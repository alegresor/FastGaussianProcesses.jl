var documenterSearchIndex = {"docs":
[{"location":"#Fast-Gaussian-Processes","page":"Home","title":"Fast Gaussian Processes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"See the Tutorial for instructions on how to use this package. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"FastGaussianProcesses.jl implements fast construction methods for Gaussian Process regression when one has control over the design of experiments. By matching shift invariant kernels to N quasi-random sampling locations, the cost of GP fitting (including optimization of kernel parameters) is reduced from  mathcalO(N^3) to mathcalO(N log N). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"FastGaussianProcesses.jl also includes methods which quickly incorporate gradient information into the GP model. Suppose we want to fit a GP to f 01^s to mathbbR and we have access to f^(boldsymbolbeta_j) for 1 leq j leq M, perhaps via automatic differentiation. For example, if we have access to f and nabla f then M=s+1 and the derivative orders are boldsymbolbeta_s+1 = boldsymbol0 and boldsymbolbeta_j = boldsymbole_j for 1 leq j leq s where boldsymbole_j is 1 at index j and 0 everywhere else. In the generalized case, the cost of GP fitting, including optimization of kernel parameters, is reduced from  mathcalO(M^3N^3) to mathcalO(M^2 N log N + M^3 N). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package implements three flavors of GPs which we discuss in the next three sections. \"Fast\" GP indicates this is a reduced cost method requiring control over the design of experiments and using sampling sequences in 01^s. QMCGenerators.jl is used to generate these lattice and digital sequences in base 2. The GPs are compared on the example function f01 to mathbbR defined as  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"f(x) = f^(0)(x) = xsin(10x)","category":"page"},{"location":"","page":"Home","title":"Home","text":"with gradient ","category":"page"},{"location":"","page":"Home","title":"Home","text":"f^(1)(x) = sin(10x)+10xcos(10x)","category":"page"},{"location":"#Fast-GP-Lattice","page":"Home","title":"Fast GP Lattice","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fast GP with a lattice sampling sequence (base 2) and matching shift invariant kernel. This additionally supports derivative information beyond the gradient i.e. second derivatives and beyond. Kernels of varying derivative order are shown below. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"","page":"Home","title":"Home","text":"The fitted GP is visualized its posterior mean (the blue solid line) and 95% confidence interval around each point (the shaded light blue region). Notice the periodicity implied by the above kernels is assumed by the GP.  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"#Fast-GP-Digital","page":"Home","title":"Fast GP Digital","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Fast GP with a digital sampling sequence (base 2) and matching digitally shift invariant kernel. Notice the discontinuities in the kernels.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"#GP-RBF","page":"Home","title":"GP RBF","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Slow) GP with arbitrary samples in mathbbR^s and a radial basis function (RBF) kernel. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: image)","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Jagadeeswaran, R. (2019). Fast automatic Bayesian cubature using matching kernels and designs. Illinois Institute of Technology.\nJagadeeswaran, R., & Hickernell, F. J. (2019). Fast automatic Bayesian cubature using lattice sampling. Statistics and Computing, 29(6), 1215-1229.\nJagadeeswaran, R., & Hickernell, F. J. (2022). Fast Automatic Bayesian Cubature Using Sobol’ Sampling. In Advances in Modeling and Simulation: Festschrift for Pierre L'Ecuyer (pp. 301-318). Cham: Springer International Publishing.\nKaarnioja, V., Kuo, F. Y., & Sloan, I. H. (2023). Lattice-based kernel approximation and serendipitous weights for parametric PDEs in very high dimensions. arXiv preprint arXiv:2303.17755.\nRasmussen, C. E., & Williams, C. K. (2006). Gaussian processes for machine learning (Vol. 1, p. 159). Cambridge, MA: MIT press.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"DocTestSetup = quote\n    using FastGaussianProcesses\n    using QMCGenerators\n    using Printf\n    using CairoMakie\n    Base.show(io::IO, f::Float64) = @printf(io, \"%1.3f\", f)\nend","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Pages = [\"tutorial.md\"]\nDepth = 4","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To begin, install this package with ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Pkg.add(\"FastGaussianProcesses\")","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"and then import via","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using FastGaussianProcesses","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To construct specific sampling sequences we also need to be ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using QMCGenerators","category":"page"},{"location":"tutorial/#Kernel-1s-Lines","page":"Tutorial","title":"Kernel 1s Lines","text":"","category":"section"},{"location":"tutorial/#Lattice-Kernels","page":"Tutorial","title":"Lattice Kernels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"PLOTDIR = joinpath(@__DIR__,\"src/assets\")\nplot_gp_kernel_latticeseqb2_1s_lines(figpath=joinpath(PLOTDIR,\"kernels.lines.latticeseqb2.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Digital-Kernels","page":"Tutorial","title":"Digital Kernels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot_gp_kernel_digitalseqb2g_1s_lines(figpath=joinpath(PLOTDIR,\"kernels.lines.digitalseqb2g.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#RBF-Kernels","page":"Tutorial","title":"RBF Kernels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot_gp_kernel_rbf_1s_lines(figpath=joinpath(PLOTDIR,\"kernels.lines.rbf.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Kernel-1s-Contours-and-Surfaces","page":"Tutorial","title":"Kernel 1s Contours and Surfaces","text":"","category":"section"},{"location":"tutorial/#Lattice-Kernels-2","page":"Tutorial","title":"Lattice Kernels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"PLOTDIR = joinpath(@__DIR__,\"src/assets\")\nplot_gp_kernel_latticeseqb2_1s_contsurfs(figpath=joinpath(PLOTDIR,\"kernels.contsurf.latticeseqb2.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Digital-Kernels-2","page":"Tutorial","title":"Digital Kernels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot_gp_kernel_digitalseqb2g_1s_contsurfs(figpath=joinpath(PLOTDIR,\"kernels.contsurf.digitalseqb2g.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#RBF-Kernels-2","page":"Tutorial","title":"RBF Kernels","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot_gp_kernel_rbf_1s_contsurfs(figpath=joinpath(PLOTDIR,\"kernels.contsurf.rbf.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#GPs-1s","page":"Tutorial","title":"GPs 1s","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"function f1s(x::Vector{Float64})\n    [   x[1]*sin(10*x[1]), # f^{(0)}\n        sin(10*x[1]) + 10*x[1]*cos(10*x[1]) # f^{(1)}\n    ]\nend","category":"page"},{"location":"tutorial/#GP-Lattice-Kernel","page":"Tutorial","title":"GP Lattice Kernel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gp_1s_latticeseqb2 = FastGaussianProcess(f1s,RandomShift(LatticeSeqB2(1),1,7),2^2;β=[i for i=0:1,j=1:1],verbose=0)\nplot_gp_optimization(gp_1s_latticeseqb2,figpath=joinpath(PLOTDIR,\"optim.1s.latticeseqb2.svg\"))\nplot_gp_1s(gp_1s_latticeseqb2;f=f1s,β=[0,1],figpath=joinpath(PLOTDIR,\"gp.1s.latticeseqb2.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#GP-Digital-Kernel","page":"Tutorial","title":"GP Digital Kernel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gp_1s_digitalseqb2g = FastGaussianProcess(f1s,RandomDigitalShift(DigitalSeqB2G(1),1,7),2^2;β=[i for i=0:1,j=1:1],verbose=0)\nplot_gp_optimization(gp_1s_digitalseqb2g,figpath=joinpath(PLOTDIR,\"optim.1s.digitalseqb2g.svg\"))\nplot_gp_1s(gp_1s_digitalseqb2g;f=f1s,β=[0,1],figpath=joinpath(PLOTDIR,\"gp.1s.digitalseqb2g.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#GP-RBF-Kernel","page":"Tutorial","title":"GP RBF Kernel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gp_1s_rbf = GaussianProcessRBF(f1s,IIDU01Seq(1,7),2^2;β=[i for i=0:1,j=1:1],verbose=0)\nplot_gp_optimization(gp_1s_rbf,figpath=joinpath(PLOTDIR,\"optim.1s.rbf.svg\"))\nplot_gp_1s(gp_1s_rbf;f=f1s,β=[0,1],figpath=joinpath(PLOTDIR,\"gp.1s.rbf.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#GPs-2s","page":"Tutorial","title":"GPs 2s","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"function f2s(x::Vector{Float64})\n    [   x[1]*sin(10*x[1]) - x[2]*cos(10*x[2]), # f^{(0,0)},\n        sin(10*x[1]) + 10*x[1]*cos(10*x[1]), # f^{(1,0)}\n        10*x[2]*sin(10*x[2]) - cos(10*x[2]) # f^{(0,1)}\n    ]\nend","category":"page"},{"location":"tutorial/#GP-Lattice-Kernel-2","page":"Tutorial","title":"GP Lattice Kernel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gp_2s_latticeseqb2 = FastGaussianProcess(f2s,RandomShift(LatticeSeqB2(2),1,7),2^8;β=[0 0; 1 0; 0 1],optim_steps=40,verbose=0)\nplot_gp_optimization(gp_2s_latticeseqb2,figpath=joinpath(PLOTDIR,\"optim.2s.latticeseqb2.svg\"))\nplot_gp_2s(gp_2s_latticeseqb2;f=f2s,β=[0 0; 1 0; 0 1],figpath=joinpath(PLOTDIR,\"gp.2s.latticeseqb2.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#GP-Digital-Kernel-2","page":"Tutorial","title":"GP Digital Kernel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gp_2s_digitalseqb2g = FastGaussianProcess(f2s,RandomDigitalShift(DigitalSeqB2G(2),1,7),2^8;β=[0 0; 1 0; 0 1],optim_steps=40,verbose=0)\nplot_gp_optimization(gp_2s_digitalseqb2g,figpath=joinpath(PLOTDIR,\"optim.2s.digitalseqb2g.svg\"))\nplot_gp_2s(gp_2s_digitalseqb2g;f=f2s,β=[0 0; 1 0; 0 1],figpath=joinpath(PLOTDIR,\"gp.2s.digitalseqb2g.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#GP-RBF-Kernel-2","page":"Tutorial","title":"GP RBF Kernel","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gp_2s_rbf = GaussianProcessRBF(f2s,IIDU01Seq(2,7),2^8;β=[0 0; 1 0; 0 1],optim_steps=40,verbose=0)\nplot_gp_optimization(gp_2s_rbf,figpath=joinpath(PLOTDIR,\"optim.2s.rbf.svg\"))\nplot_gp_2s(gp_2s_rbf;f=f2s,β=[0 0; 1 0; 0 1],figpath=joinpath(PLOTDIR,\"gp.2s.rbf.svg\"))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"},{"location":"tutorial/#Logo","page":"Tutorial","title":"Logo","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In the following we must be ","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using CairoMakie","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"m = 2; n = 2^m\nx = FirstLinear(RandomDigitalShift(DigitalSeqB2G(1),1,11),m)[:,1]\nkmat = Matrix{Float64}(undef,2*n,2*n)\nfor i1=1:n,i2=1:n,β1=0:1,β2=0:1 kmat[β1*n+i1,β2*n+i2] = ((β1+β2)==0)+kernel_digshiftinvar_s1(x[i1],x[i2],β1,β2,4) end \nfig = Figure(resolution=(400,400),backgroundcolor=:transparent)\nax = Axis(fig[1,1],backgroundcolor=:transparent)\nheatmap!(ax,1:2*n,2*n:-1:1,kmat,colormap=:julia_colorscheme)#:nipy_spectral) # https://docs.juliahub.com/MakieGallery/Ql23q/0.2.17/generated/colors.html\nhidespines!(ax); hidedecorations!(ax); hidexdecorations!(ax,grid = false); hideydecorations!(ax, ticks = false)\nsave(joinpath(PLOTDIR,\"logo.svg\"),fig)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: image)","category":"page"}]
}
