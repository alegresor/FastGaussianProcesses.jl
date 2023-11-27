module FastGaussianProcesses

using QMCGenerators
import FFTW: fft,ifft
import Hadamard: fwht_natural,ifwht_natural
import LinearAlgebra: logdet,tr
import Printf: @printf
import CairoMakie
using LaTeXStrings
import Distributions: quantile,Normal

include("util.jl")
export RationalToBinary,Float64ToBinary

include("gp.latticeseqb2.jl")
export kernel_shiftinvar,GaussianProcessLatticeSeqB2

include("gp.digitalseqb2g.jl")
export kernel_digshiftinvar

export mean_post,cov_post,var_post

include("plots.jl")
export plot_gp_kernel_latticeseqb2_1s_lines,plot_gp_kernel_digitalseqb2g_1s_lines,plot_gp_kernel_shiftinvar_1s_contour,plot_gp_optimization,plot_gp_1s,plot_gp_2s

end