function k4sumterm(xb::UInt64,t::Int64)
    total = BigFloat(0.) 
    for a1=0:t-1
        total += (-1)^((xb>>(t-a1-1))&1)/BigFloat(2)^(3*a1)
    end 
    Float64(total)
end 

KTILDE = Dict(
    2 => (β::Int64,x::Float64,xb::UInt64,t::Int64) -> x==0 ? 5/2 : -β*x + 5/2*(1-2. ^(-β)),
    3 => (β::Int64,x::Float64,xb::UInt64,t::Int64) -> x==0 ? 43/18 : β*x^2 - 5*(1-2. ^(-β))*x + 43/18*(1-2. ^(-2*β)),
    4 => (β::Int64,x::Float64,xb::UInt64,t::Int64) -> x==0 ? 701/294 : -2/3*β*x^3 + 5*(1-2. ^(-β))*x^2 - 43/9*(1-2. ^(-2*β))*x + 701/294*(1-2. ^(-3*β)) + β*(1/48*k4sumterm(xb,t)-1/42))

function kernel_digshiftinvar_s1(x1b::UInt64,x2b::UInt64,β1::Int64,β2::Int64,α::Int64,t::Int64)
    xb = x1b ⊻ x2b
    x = xb * 2. ^(-t)
    β = x == 0 ? 0 : Int64(-floor(log2(x)))
    ktildeidx = α-β1-β2
    (-2)^(β1+β2)*KTILDE[ktildeidx](β,x,xb,t)
end 

function kernel_digshiftinvar(x1b::Vector{UInt64},x2b::Vector{UInt64},β1::Vector{Int64},β2::Vector{Int64},α::Int64,γ::Float64,η::Vector{Float64},s::Int64,t::Int64)
    kvals = map(j->kernel_digshiftinvar_s1(x1b[j],x2b[j],β1[j],β2[j],α,t),1:s)
    γ*prod(map(j->β1[j]+β2[j] == 0 ? 1+η[j]*kvals[j] : η[j]*kvals[j],1:s))
end
kernel_digshiftinvar(x1b::Vector{UInt64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},α::Int64,γ::Float64,η::Vector{Float64},s::Int64,t::Int64) = kernel_digshiftinvar(x1b,Float64ToBinary.(x2,t),β1,β2,α,γ,η,s,t)
kernel_digshiftinvar(x1::Vector{Float64},x2b::Vector{UInt64},β1::Vector{Int64},β2::Vector{Int64},α::Int64,γ::Float64,η::Vector{Float64},s::Int64,t::Int64) = kernel_digshiftinvar(Float64ToBinary.(x1,t),x2b,β1,β2,α,γ,η,s,t)
kernel_digshiftinvar(x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},α::Int64,γ::Float64,η::Vector{Float64},s::Int64,t::Int64) = kernel_digshiftinvar(Float64ToBinary.(x1,t),Float64ToBinary.(x2,t),β1,β2,α,γ,η,s,t)
kernel_digshiftinvar(x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},α::Int64,γ::Float64,η::Vector{Float64},s::Int64) = kernel_digshiftinvar(x1,x2,β1,β2,α,γ,η,s,53)

mutable struct GaussianProcessDigitalSeqB2G
    s::Int64
    n::Int64
    α::Int64
    γ::Float64
    η::Vector{Float64}
    β::Matrix{Int64}
    n_β::Int64
    xb::Matrix{UInt64}
    x::Matrix{Float64}
    y::Matrix{Float64}
    ytilde::Matrix{Float64}
    ktilde::Array{Float64}
    coeffs::Matrix{Float64}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
    t::Int64
end

function GaussianProcessDigitalSeqB2G(f::Function,ds::Union{DigitalSeqB2G,RandomDigitalShift},n::Int64;β::Union{Nothing,Matrix{Int64}}=nothing,α::Int64=4,γ::Float64=1.,η::Union{Float64,Vector{Float64}}=1.,ζ::Union{Float64,Vector{Float64}}=1e-16,optim_steps::Int64=100,learningrate::Float64=1e-1,decayrate::Float64=.9,verbose::Int64=10)
    verbosebool = verbose > 0
    @assert log2(n)%1==0
    m = Int64(log2(n))
    t = ds.t
    xb = FirstLinearBinary(ds,m)
    x = BinaryToFloat64(xb,ds)
    s = size(x,2)
    if β===nothing β = zeros(Int64,1,s) end 
    n_β = size(β,1)
    @assert size(β,2)==s
    @assert 2*α in keys(BERNOULLIPOLYS)
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,n_β) else @assert size(ζ)==(n_β,) end 
    p = 1+s+n_β
    y = reshape(vcat([f(x[i,:]) for i=1:n]'...),n,n_β)
    ytilde = sqrt(n).*fwht_natural(y,1)
    losses,γs,ηs,ζs = zeros(Float64,optim_steps+1),zeros(Float64,optim_steps+1),zeros(Float64,optim_steps+1,s),zeros(Float64,optim_steps+1,n_β)
    if verbosebool println("QGP Optimization Loss") end
    Δ = zeros(Float64,p)
    for step=1:optim_steps
        γs[step] = γ; ηs[step,:] = η; ζs[step,:] = ζ
        kvals = [kernel_digshiftinvar(xb[i,:],xb[1,:],β[k,:],β[l,:],α,γ,η,s,t) for i=1:n,k=1:n_β,l=1:n_β]
        kvalsnoisy = copy(kvals); for k=1:n_β kvalsnoisy[1,k,k] = kvalsnoisy[1,k,k]+ζ[k] end
        ktilde = n.*fwht_natural(kvalsnoisy,1)
        coeffs_perm = permutedims(hcat(map(i->ktilde[i,:,:]\ytilde[i,:],1:n)...))
        losses[step] = sum(map(i->logdet(ktilde[i,:,:]),1:n))+sum(conj.(ytilde).*coeffs_perm)
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,losses[step]) end
        ∂k∂γ = kvals./γ
        ∂k∂η = zeros(Float64,s,n,n_β,n_β)
        for j=1:s,i=1:n,k=1:n_β,l=1:n_β
            po_kj,po_lj = β[k,j],β[l,j]
            kval_jikl = kernel_digshiftinvar_s1(xb[i,j],xb[1,j],po_kj,po_lj,α,t)
            if kval_jikl==0 ∂k∂η[j,i,k,l] = 0; continue end 
            denom = po_kj+po_lj == 0 ? 1+η[j]*kval_jikl : η[j]*kval_jikl
            ∂k∂η[j,i,k,l] = kvals[i,k,l]*kval_jikl/denom
        end
        ∂k∂ζ = zeros(Float64,n_β,n,n_β,n_β); for k=1:n_β ∂k∂ζ[k,1,k,k] = 1 end
        θ = [γ,η...,ζ...]
        ∂k∂θ = zeros(Float64,p,n,n_β,n_β); ∂k∂θ[1,:,:,:] = ∂k∂γ; ∂k∂θ[2:s+1,:,:,:] = ∂k∂η; ∂k∂θ[s+2:end,:,:,:] = ∂k∂ζ
        ∂ktilde∂θ = n.*fwht_natural(∂k∂θ,2)
        ∂L∂θ = [sum(map(i->tr(ktilde[i,:,:]\∂ktilde∂θ[j,i,:,:])-coeffs_perm[i,:]'*∂ktilde∂θ[j,i,:,:]*coeffs_perm[i,:],1:n)) for j=1:p]
        ∂L∂logθ = ∂L∂θ.*θ # since ∂θ∂logθ = θ
        Δ = decayrate*Δ.+(1-decayrate)*∂L∂logθ.^2
        θ = exp.(log.(θ)-learningrate*Δ.^(-1/2).*∂L∂logθ)
        γ = θ[1]; η = θ[2:1+s]; ζ = θ[2+s:end]
    end
    ktilde = sqrt(n).*fwht_natural([kernel_digshiftinvar(xb[i,:],xb[1,:],β[k,:],β[l,:],α,γ,η,s,t)+ζ[k]*(k==l)*(i==1) for i=1:n,k=1:n_β,l=1:n_β],1)
    coeffs_perm = permutedims(hcat(map(i->ktilde[i,:,:]\ytilde[i,:],1:n)...))
    losses[end] = sum(map(i->logdet(ktilde[i,:,:]),1:n))+sum(conj.(ytilde).*coeffs_perm)
    γs[end] = γ; ηs[end,:] = η; ζs[end,:] = ζ
    coeffs = ifwht_natural(coeffs_perm,1)/n
    GaussianProcessDigitalSeqB2G(s,n,α,γ,η,β,n_β,xb,x,y,ytilde,ktilde,coeffs,losses,γs,ηs,ζs,t)
end
GaussianProcessDigitalSeqB2G(f::Function,s::Int64,n::Int64;kwargs...) = GaussianProcessDigitalSeqB2G(f,RandomDigitalShift(DigitalSeqB2G(LinearMatrixScramble(s))),n;kwargs...)

function mean_post(gp::GaussianProcessDigitalSeqB2G,x::Union{Vector{UInt64},Vector{Float64}},β::Vector{Int64})
    kmat = [kernel_digshiftinvar(x,gp.xb[i,:],β,gp.β[j,:],gp.α,gp.γ,gp.η,gp.s,gp.t) for i=1:gp.n,j=1:gp.n_β]
    sum(gp.coeffs.*kmat)
end 
(gp::GaussianProcessDigitalSeqB2G)(x::Union{Vector{UInt64},Vector{Float64}},β::Vector{Int64}) = mean_post(gp,x,β)

function cov_post(gp::GaussianProcessDigitalSeqB2G,x1::Union{Vector{UInt64},Vector{Float64}},x2::Union{Vector{UInt64},Vector{Float64}},β1::Vector{Int64},β2::Vector{Int64})
    kval = kernel_digshiftinvar(x1,x2,β1,β2,gp.α,gp.γ,gp.η,gp.s,gp.t)
    k1vec = [kernel_digshiftinvar(x1,gp.x[i,:],β1,gp.β[j,:],gp.α,gp.γ,gp.η,gp.s,gp.t) for i=1:gp.n,j=1:gp.n_β]
    k2vec = [kernel_digshiftinvar(x2,gp.x[i,:],β2,gp.β[j,:],gp.α,gp.γ,gp.η,gp.s,gp.t) for i=1:gp.n,j=1:gp.n_β]
    k2vectilde = sqrt(gp.n)*fwht_natural(k2vec,1)
    coeffs = ifwht_natural(permutedims(hcat(map(i->gp.ktilde[i,:,:]\k2vectilde[i,:],1:gp.n)...)),1)/gp.n
    kval-sum(k1vec.*coeffs)
end

function var_post(gp::GaussianProcessDigitalSeqB2G,x::Union{Vector{UInt64},Vector{Float64}},β::Vector{Int64})
    var = cov_post(gp,x,x,β,β)
    -1e-12<var<0 ? 0. : var
end