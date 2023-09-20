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
    γ*prod(map(j->partial_x1_order[j]+partial_x2_order[j] == 0 ? 1+η[j]*kvals[j] : η[j]*kvals[j],1:s))
end

mutable struct GaussianProcessLatticeSeqB2
    s::Int64
    f::Function 
    n::Int64
    o::Int64
    γ::Float64
    η::Vector{Float64}
    partial_orders::Matrix{Int64}
    r::Int64
    x::Matrix{Float64}
    y::Matrix{Float64}
    ytilde::Matrix{ComplexF64}
    ktilde::Array{ComplexF64}
    coeffs::Matrix{Float64}
end

function GaussianProcessLatticeSeqB2(f::Function,n::Int64,ls::Union{LatticeSeqB2,RandomShift},o::Int64,γ::Float64,η::Vector{Float64},ζ::Vector{Float64},partial_orders::Matrix{Int64},optim_steps::Int64,verbose::Int64)
    verbosebool = verbose > 0
    @assert log2(n)%1==0
    m = Int64(log2(n))
    x = FirstLinear(ls,m)
    s = size(x,2)
    r = size(partial_orders,1)
    p = 1+s+r
    @assert size(partial_orders,2)==s
    y = [f(x[i,:],partial_orders[l,:]) for i=1:n,l=1:r]
    ytilde = fft(y,1)
    losses,γs,ηs,ζs = zeros(Float64,optim_steps+1),zeros(Float64,optim_steps+1),zeros(Float64,optim_steps+1,s),zeros(Float64,optim_steps,r)
    if verbosebool println("QGP Optimization Loss") end 
    for step=1:optim_steps
        γs[step] = γ; ηs[step,:] = η; ζs[step,:] = ζ
        kvals = [kernel_lattice(x[i,:],x[1,:],o,γ,η,partial_orders[k,:],partial_orders[l,:],s) for i=1:n,k=1:r,l=1:r]
        kvalsnoisy = copy(kvals); for k=1:r kvalsnoisy[1,k,k] = kvalsnoisy[1,k,k]+ζ[k] end
        ktilde = fft(kvalsnoisy,1)
        coeffs_perm = permutedims(hcat(map(i->ktilde[i,:,:]\ytilde[i,:],1:n)...))
        losses[step] = real.(sum(map(i->logdet(ktilde[i,:,:]),1:n))+sum(conj.(ytilde).*coeffs_perm))
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,losses[step]) end
        ∂k∂γ = kvals./γ
        ∂k∂η = zeros(Float64,s,n,r,r)
        for j=1:s,i=1:n,k=1:r,l=1:r
            po_kj,po_lj = partial_orders[k,j],partial_orders[l,j]
            kval_jikl = kernel_lattice_s1(x[i,j],x[1,j],po_kj,po_lj,o)
            if kval_jikl==0 ∂k∂η[j,i,k,l] = 0; continue end 
            denom = po_kj+po_lj == 0 ? 1+η[j]*kval_jikl : η[j]*kval_jikl
            ∂k∂η[j,i,k,l] = kvals[i,k,l]*kval_jikl/denom
        end
        ∂k∂ζ = zeros(Float64,r,n,r,r); for k=1:r ∂k∂ζ[k,1,k,k] = 1 end
        θ = [γ,η...,ζ...]
        ∂k∂θ = zeros(Float64,p,n,r,r); ∂k∂θ[1,:,:,:] = ∂k∂γ; ∂k∂θ[2:s+1,:,:,:] = ∂k∂η; ∂k∂θ[s+2:end,:,:,:] = ∂k∂ζ
        ∂ktilde∂θ = fft(∂k∂θ,2)
        ∂L∂θ = real.([sum(map(i->tr(ktilde[i,:,:]\∂ktilde∂θ[j,i,:,:])-coeffs_perm[i,:]'*∂ktilde∂θ[j,i,:,:]*coeffs_perm[i,:],1:n)) for j=1:p])
        ∂L∂logθ = ∂L∂θ.*θ # since ∂θ∂logθ = θ
        θ = exp.(log.(θ) - .05 .* ∂L∂logθ) # TODO
        γ = θ[1]; η = θ[2:1+s]; ζ = θ[2+s:end]
    end
    ktilde = fft([kernel_lattice(x[i,:],x[1,:],o,γ,η,partial_orders[k,:],partial_orders[l,:],s)+ζ[k]*(k==l) for i=1:n,k=1:r,l=1:r],1)
    coeffs = real.(ifft(permutedims(hcat(map(i->ktilde[i,:,:]\ytilde[i,:],1:n)...)),1))
    GaussianProcessLatticeSeqB2(s,f,n,o,γ,η,partial_orders,r,x,y,ytilde,ktilde,coeffs)
end

function mean_post(gp::GaussianProcessLatticeSeqB2,x::Vector{Float64},partial_order::Vector{Int64})
    kmat = [kernel_lattice(x,gp.x[i,:],gp.o,gp.γ,gp.η,partial_order,gp.partial_orders[j,:],gp.s) for i=1:gp.n,j=1:gp.r]
    sum(gp.coeffs.*kmat)
end 

(gp::GaussianProcessLatticeSeqB2)(x::Vector{Float64},partial_order::Vector{Int64}) = mean_post(gp,x,partial_order)

function cov_post(gp::GaussianProcessLatticeSeqB2,x1::Vector{Float64},x2::Vector{Float64},partial_x1_order::Vector{Int64},partial_x2_order::Vector{Int64})
    kval = kernel_lattice(x1,x2,gp.o,gp.γ,gp.η,partial_x1_order,partial_x2_order,gp.s)
    k1vec = [kernel_lattice(x1,gp.x[i,:],gp.o,gp.γ,gp.η,partial_x1_order,gp.partial_orders[j,:],gp.s) for i=1:gp.n,j=1:gp.r]
    k2vec = [kernel_lattice(x2,gp.x[i,:],gp.o,gp.γ,gp.η,partial_x2_order,gp.partial_orders[j,:],gp.s) for i=1:gp.n,j=1:gp.r]
    k2vectilde = fft(k2vec,1)
    coeffs = real.(ifft(permutedims(hcat(map(i->gp.ktilde[i,:,:]\k2vectilde[i,:],1:gp.n)...)),1))
    kval-sum(k1vec.*coeffs)
end

var_post(gp::GaussianProcessLatticeSeqB2,x::Vector{Float64},partial_order::Vector{Int64}) = cov_post(gp,x,x,partial_order,partial_order)