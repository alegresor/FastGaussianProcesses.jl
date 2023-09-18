module FastGaussianProcesses

using QMCGenerators
import FFTW: fft, ifft
import Hadamard: fwht_natural,ifwht_natural

include("gp.latticeseqb2.jl")
export kernel_lattice,GaussianProcessLatticeSeqB2

export mean_post

end