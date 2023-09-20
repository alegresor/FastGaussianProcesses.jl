module FastGaussianProcesses

using QMCGenerators
import FFTW: fft, ifft
import Hadamard: fwht_natural,ifwht_natural
import LinearAlgebra: logdet,tr
import Printf: @printf

include("gp.latticeseqb2.jl")
export kernel_lattice,GaussianProcessLatticeSeqB2

export mean_post,cov_post,var_post

end