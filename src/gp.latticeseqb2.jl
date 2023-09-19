BERNOULLIPOLYS = Dict(
    0 =>  x -> 1,
    1 =>  x -> x-1/2,
    2 =>  x -> x^2-x+1/6,
    3 =>  x -> x^3-3/2*x^2+1/2*x,
    4 =>  x -> x^4-2*x^3+x^2-1/30,
    5 =>  x -> x^5-5/2*x^4+5/3*x^3-1/6*x,
    6 =>  x -> x^6-3*x^5+5/2*x^4-1/2*x^2+1/42,
    7 =>  x -> x^7-7/2*x^6+7/2*x^5-7/6*x^3+1/6*x,
    8 =>  x -> x^8-4*x^7+14/3*x^6-7/3*x^4+2/3*x^2-1/30,
    9 =>  x -> x^9-9/2*x^8+6*x^7-21/5*x^5+2*x^3-3/10*x,
    10 => x -> x^10-5*x^9+15/2*x^8-7*x^6+5*x^4-3/2*x^2+5/66)

function kernel_lattice_s1(x1::Float64,x2::Float64,partial_x1_order::Int64,partial_x2_order::Int64,o::Int64)
    bpolyidx = 2*o-partial_x1_order-partial_x2_order
    (-1)^(o+1+partial_x2_order)*prod((bpolyidx+1):(2*o))*BERNOULLIPOLYS[bpolyidx](mod(x1-x2,1))
end 

function kernel_lattice(x1::Vector{Float64},x2::Vector{Float64},o::Int64,γ::Float64,η::Vector{Float64},partial_x1_order::Vector{Int64},partial_x2_order::Vector{Int64},s::Int64)
    kvals = map(j->kernel_lattice_s1(x1[j],x2[j],partial_x1_order[j],partial_x2_order[j],o),1:s)
    γ*prod(map(j-> partial_x1_order[j]+partial_x2_order[j] == 0 ? 1+η[j]*kvals[j] : η[j]*kvals[j],1:s))
end

mutable struct GaussianProcessLatticeSeqB2
    s::Int64
    f::Function 
    n::Int64
    o::Int64
    γ::Float64
    η::Vector{Float64}
    partial_orders::Matrix{Int64}
    n_orders::Int64
    x::Matrix{Float64}
    y::Matrix{Float64}
    ytilde::Matrix{ComplexF64}
    ktilde::Array{ComplexF64}
    coeffs::Matrix{Float64}
end

function GaussianProcessLatticeSeqB2(f::Function,n::Int64,ls::Union{LatticeSeqB2,RandomShift},o::Int64,γ::Float64,η::Vector{Float64},partial_orders::Matrix{Int64})
    @assert log2(n)%1==0
    m = Int64(log2(n))
    x = FirstLinear(ls,m)
    s = size(x,2)
    n_orders = size(partial_orders,1)
    @assert size(partial_orders,2)==s
    y = [f(x[i,:],partial_orders[j,:]) for i=1:n,j=1:n_orders]
    ytilde = fft(y,1)
    ktilde = zeros(ComplexF64,n,n_orders,n_orders)
    for k=0:n_orders-1
        for j=0:n_orders-1 # TODO use symmetry and real only / imag only to make more efficient
            ktilde[:,k+1,j+1] = fft([kernel_lattice(x[i,:],x[1,:],o,γ,η,partial_orders[k+1,:],partial_orders[j+1,:],s) for i=1:n])
        end
    end
    freqs = zeros(ComplexF64,n,n_orders)
    for i=1:n freqs[i,:] = ktilde[i,:,:]\ytilde[i,:] end 
    coeffs = real.(ifft(freqs,1))
    GaussianProcessLatticeSeqB2(s,f,n,o,γ,η,partial_orders,n_orders,x,y,ytilde,ktilde,coeffs)
end

function mean_post(gp::GaussianProcessLatticeSeqB2,x::Vector{Float64},partial_order::Vector{Int64})
    kmat = [kernel_lattice(x,gp.x[i,:],gp.o,gp.γ,gp.η,partial_order,gp.partial_orders[j,:],gp.s) for i=1:gp.n,j=1:gp.n_orders]
    sum(gp.coeffs.*kmat)
end 

(gp::GaussianProcessLatticeSeqB2)(x::Vector{Float64},partial_order::Vector{Int64}) = mean_post(gp,x,partial_order)

function cov_post(gp::GaussianProcessLatticeSeqB2,x1::Vector{Float64},x2::Vector{Float64},partial_x1_order::Vector{Int64},partial_x2_order::Vector{Int64})
    kval = kernel_lattice(x1,x2,gp.o,gp.γ,gp.η,partial_x1_order,partial_x2_order,gp.s)
    k1vec = [kernel_lattice(x1,gp.x[i,:],gp.o,gp.γ,gp.η,partial_x1_order,gp.partial_orders[j,:],gp.s) for i=1:gp.n,j=1:gp.n_orders]
    k2vec = [kernel_lattice(x2,gp.x[i,:],gp.o,gp.γ,gp.η,partial_x2_order,gp.partial_orders[j,:],gp.s) for i=1:gp.n,j=1:gp.n_orders]
    k2vectilde = fft(k2vec,1)
    freqs = zeros(ComplexF64,gp.n,gp.n_orders)
    for i=1:gp.n freqs[i,:] = gp.ktilde[i,:,:]\k2vectilde[i,:] end 
    coeffs = real.(ifft(freqs,1))
    kval-sum(k1vec.*coeffs)
end

var_post(gp::GaussianProcessLatticeSeqB2,x::Vector{Float64},partial_order::Vector{Int64}) = cov_post(gp,x,x,partial_order,partial_order)