module FastGaussianProcesses

using QMCGenerators
import FFTW: fft
import Hadamard: fwht_natural
import LinearAlgebra: logdet,tr
import Printf: @printf
import CairoMakie
using LaTeXStrings
import Distributions: quantile,Normal

include("util.jl")
export RationalToBinary,Float64ToBinary

include("gp.jl")
export FastGaussianProcess,mean_post,cov_post,var_post

include("gp.latticeseqb2.jl")
export kernel_shiftinvar,GaussianProcessLatticeSeqB2

include("gp.digitalseqb2g.jl")
export kernel_digshiftinvar,GaussianProcessDigitalSeqB2G

include("plots.jl")
export plot_gp_kernel_latticeseqb2_1s_lines,plot_gp_kernel_digitalseqb2g_1s_lines,plot_gp_kernel_latticeseqb2_1s_contsurfs,plot_gp_kernel_digitalseqb2g_1s_contsurfs,plot_gp_optimization,plot_gp_1s,plot_gp_2s

end