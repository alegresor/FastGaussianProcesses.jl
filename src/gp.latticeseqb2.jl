BERNOULLIPOLYS = Dict(
    #0 =>  x -> 1,
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

function kernel_lattice_s1(x1::Float64,x2::Float64,β1::Int64,β2::Int64,α::Int64)
    bpolyidx = 2*α-β1-β2
    (-1)^(α+1+β2)*prod((bpolyidx+1):(2*α))*BERNOULLIPOLYS[bpolyidx](mod(x1-x2,1))
end 

function kernel_lattice(x1::Vector{Float64},x2::Vector{Float64},α::Int64,γ::Float64,η::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},s::Int64)
    kvals = map(j->kernel_lattice_s1(x1[j],x2[j],β1[j],β2[j],α),1:s)
    γ*prod(map(j->β1[j]+β2[j] == 0 ? 1+η[j]*kvals[j] : η[j]*kvals[j],1:s))
end

mutable struct GaussianProcessLatticeSeqB2
    s::Int64
    f::Function 
    n::Int64
    α::Int64
    γ::Float64
    η::Vector{Float64}
    β::Matrix{Int64}
    r::Int64
    x::Matrix{Float64}
    y::Matrix{Float64}
    ytilde::Matrix{ComplexF64}
    ktilde::Array{ComplexF64}
    coeffs::Matrix{Float64}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
end

function GaussianProcessLatticeSeqB2(f::Function,ls::Union{LatticeSeqB2,RandomShift},n::Int64;β::Union{Nothing,Matrix{Int64}}=nothing,α::Int64=4,γ::Float64=1.,η::Union{Float64,Vector{Float64}}=1.,ζ::Union{Float64,Vector{Float64}}=1e-16,optim_steps::Int64=100,learningrate::Float64=1e-1,decayrate::Float64=.9,verbose::Int64=10)
    verbosebool = verbose > 0
    @assert log2(n)%1==0
    m = Int64(log2(n))
    x = FirstLinear(ls,m)
    s = size(x,2)
    if β===nothing β = zeros(Int64,1,s) end 
    r = size(β,1)
    @assert size(β,2)==s
    @assert 2*α in keys(BERNOULLIPOLYS)
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,r) else @assert size(ζ)==(r,) end 
    p = 1+s+r
    y = reshape(vcat([f(x[i,:]) for i=1:n]'...),n,r)
    ytilde = fft(y,1)
    losses,γs,ηs,ζs = zeros(Float64,optim_steps+1),zeros(Float64,optim_steps+1),zeros(Float64,optim_steps+1,s),zeros(Float64,optim_steps+1,r)
    if verbosebool println("QGP Optimization Loss") end
    Δ = zeros(Float64,p)
    for step=1:optim_steps
        γs[step] = γ; ηs[step,:] = η; ζs[step,:] = ζ
        kvals = [kernel_lattice(x[i,:],x[1,:],α,γ,η,β[k,:],β[l,:],s) for i=1:n,k=1:r,l=1:r]
        kvalsnoisy = copy(kvals); for k=1:r kvalsnoisy[1,k,k] = kvalsnoisy[1,k,k]+ζ[k] end
        ktilde = fft(kvalsnoisy,1)
        coeffs_perm = permutedims(hcat(map(i->ktilde[i,:,:]\ytilde[i,:],1:n)...))
        losses[step] = real.(sum(map(i->logdet(ktilde[i,:,:]),1:n))+sum(conj.(ytilde).*coeffs_perm))
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,losses[step]) end
        ∂k∂γ = kvals./γ
        ∂k∂η = zeros(Float64,s,n,r,r)
        for j=1:s,i=1:n,k=1:r,l=1:r
            po_kj,po_lj = β[k,j],β[l,j]
            kval_jikl = kernel_lattice_s1(x[i,j],x[1,j],po_kj,po_lj,α)
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
        Δ = decayrate*Δ.+(1-decayrate)*∂L∂logθ.^2
        θ = exp.(log.(θ)-learningrate*Δ.^(-1/2).*∂L∂logθ)
        γ = θ[1]; η = θ[2:1+s]; ζ = θ[2+s:end]
    end
    ktilde = fft([kernel_lattice(x[i,:],x[1,:],α,γ,η,β[k,:],β[l,:],s)+ζ[k]*(k==l)*(i==1) for i=1:n,k=1:r,l=1:r],1)
    coeffs_perm = permutedims(hcat(map(i->ktilde[i,:,:]\ytilde[i,:],1:n)...))
    losses[end] = real.(sum(map(i->logdet(ktilde[i,:,:]),1:n))+sum(conj.(ytilde).*coeffs_perm))
    γs[end] = γ; ηs[end,:] = η; ζs[end,:] = ζ
    coeffs = real.(ifft(coeffs_perm,1))
    GaussianProcessLatticeSeqB2(s,f,n,α,γ,η,β,r,x,y,ytilde,ktilde,coeffs,losses,γs,ηs,ζs)
end

GaussianProcessLatticeSeqB2(f::Function,s::Int64,n::Int64;kwargs...) = GaussianProcessLatticeSeqB2(f,RandomShift(LatticeSeqB2(s)),n;kwargs...)

function mean_post(gp::GaussianProcessLatticeSeqB2,x::Vector{Float64},β::Vector{Int64})
    kmat = [kernel_lattice(x,gp.x[i,:],gp.α,gp.γ,gp.η,β,gp.β[j,:],gp.s) for i=1:gp.n,j=1:gp.r]
    sum(gp.coeffs.*kmat)
end 

(gp::GaussianProcessLatticeSeqB2)(x::Vector{Float64},β::Vector{Int64}) = mean_post(gp,x,β)

function cov_post(gp::GaussianProcessLatticeSeqB2,x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64})
    kval = kernel_lattice(x1,x2,gp.α,gp.γ,gp.η,β1,β2,gp.s)
    k1vec = [kernel_lattice(x1,gp.x[i,:],gp.α,gp.γ,gp.η,β1,gp.β[j,:],gp.s) for i=1:gp.n,j=1:gp.r]
    k2vec = [kernel_lattice(x2,gp.x[i,:],gp.α,gp.γ,gp.η,β2,gp.β[j,:],gp.s) for i=1:gp.n,j=1:gp.r]
    k2vectilde = fft(k2vec,1)
    coeffs = real.(ifft(permutedims(hcat(map(i->gp.ktilde[i,:,:]\k2vectilde[i,:],1:gp.n)...)),1))
    kval-sum(k1vec.*coeffs)
end

function var_post(gp::GaussianProcessLatticeSeqB2,x::Vector{Float64},β::Vector{Int64})
    var = cov_post(gp,x,x,β,β)
    -1e-12<var<0 ? 0. : var
end