mutable struct FastGaussianProcess
    s::Int64
    n::Int64
    α::Int64
    γ::Float64
    η::Vector{Float64}
    ζ::Vector{Float64}
    β::Matrix{Int64}
    n_β::Int64
    _x::Union{Matrix{UInt64},Matrix{Float64}}
    x::Matrix{Float64}
    y::Matrix{Float64}
    λ::Union{Array{ComplexF64},Union{Float64}}
    coeffs::Matrix{Float64}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
    ft::Function # ft(a) = V^H*a
    kernel_1s::Function
end

function FastGaussianProcess(f::Function,seq::Union{LatticeSeqB2,RandomShift,DigitalSeqB2G,RandomDigitalShift},n::Int64;β::Union{Nothing,Matrix{Int64}}=nothing,α::Int64=4,γ::Float64=1.,η::Union{Float64,Vector{Float64}}=1.,ζ::Union{Float64,Vector{Float64}}=1e-16,optim_steps::Int64=100,learningrate::Float64=1e-1,decayrate::Float64=.9,verbose::Int64=10)
    @assert log2(n)%1==0
    if typeof(seq) in [LatticeSeqB2,RandomShift]
        x = _x = FirstLinear(seq,Int64(log2(n)))
        @assert 2*α in keys(BERNOULLIPOLYS)
        ft = (a,dims,n) -> conj.(fft(a,dims))/sqrt(n)
        ft_dtype = ComplexF64
        kernel_1s = kernel_shiftinvar_s1
    else # typeof(seq) in [DigitalSeqB2G,RandomDigitalShift]
        _x = FirstLinearBinary(seq,m)
        x = BinaryToFloat64(_x,seq)
        ft = (a,dims,n) -> fwht_natural(a,dims)*sqrt(n)
        ft_dtype = Float64
        kernel_1s = kernel_digshiftinvar_s1
    end 
    s = size(x,2)
    if β===nothing β = zeros(Int64,1,s) else @assert size(β,2)==s end
    n_β = size(β,1)
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,n_β) else @assert size(ζ)==(n_β,) end
    y = reshape(vcat([f(x[i,:]) for i=1:n]'...),n,n_β)
    gp = FastGaussianProcess(s,n,α,γ,η,ζ,β,n_β,_x,x,y,Array{ft_dtype}(undef,n,n_β,n_β),Matrix{Float64}(undef,n,n_β),Vector{Float64}(undef,optim_steps+1),Vector{Float64}(undef,optim_steps+1),Matrix{Float64}(undef,optim_steps+1,s),Matrix{Float64}(undef,optim_steps+1,n_β),ft,kernel_1s)
    _train(gp,verbose)
    gp 
end

function _train(gp::FastGaussianProcess,verbose::Int64)
    println("in _train")
    exit()
    verbosebool = verbose > 0
    p = 1+s+n_β
    y_ft = fft(y,1)
    if verbosebool println("QGP Optimization Loss") end
    Δ = zeros(Float64,p)
    for step=1:optim_steps
        γs[step] = γ; ηs[step,:] = η; ζs[step,:] = ζ
        kvals = [kernel_shiftinvar(x[i,:],x[1,:],β[k,:],β[l,:],α,γ,η,s) for i=1:n,k=1:n_β,l=1:n_β]
        kvalsnoisy = copy(kvals); for k=1:n_β kvalsnoisy[1,k,k] = kvalsnoisy[1,k,k]+ζ[k] end
        λ = fft(kvalsnoisy,1)
        coeffs_perm = permutedims(hcat(map(i->λ[i,:,:]\y_ft[i,:],1:n)...))
        losses[step] = real.(sum(map(i->logdet(λ[i,:,:]),1:n))+sum(conj.(y_ft).*coeffs_perm))
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,losses[step]) end
        ∂k∂γ = kvals./γ
        ∂k∂η = zeros(Float64,s,n,n_β,n_β)
        for j=1:s,i=1:n,k=1:n_β,l=1:n_β
            po_kj,po_lj = β[k,j],β[l,j]
            kval_jikl = kernel_shiftinvar_s1(x[i,j],x[1,j],po_kj,po_lj,α)
            if kval_jikl==0 ∂k∂η[j,i,k,l] = 0; continue end 
            denom = po_kj+po_lj == 0 ? 1+η[j]*kval_jikl : η[j]*kval_jikl
            ∂k∂η[j,i,k,l] = kvals[i,k,l]*kval_jikl/denom
        end
        ∂k∂ζ = zeros(Float64,n_β,n,n_β,n_β); for k=1:n_β ∂k∂ζ[k,1,k,k] = 1 end
        θ = [γ,η...,ζ...]
        ∂k∂θ = zeros(Float64,p,n,n_β,n_β); ∂k∂θ[1,:,:,:] = ∂k∂γ; ∂k∂θ[2:s+1,:,:,:] = ∂k∂η; ∂k∂θ[s+2:end,:,:,:] = ∂k∂ζ
        ∂ktilde∂θ = fft(∂k∂θ,2)
        ∂L∂θ = real.([sum(map(i->tr(λ[i,:,:]\∂ktilde∂θ[j,i,:,:])-coeffs_perm[i,:]'*∂ktilde∂θ[j,i,:,:]*coeffs_perm[i,:],1:n)) for j=1:p])
        ∂L∂logθ = ∂L∂θ.*θ # since ∂θ∂logθ = θ
        Δ = decayrate*Δ.+(1-decayrate)*∂L∂logθ.^2
        θ = exp.(log.(θ)-learningrate*Δ.^(-1/2).*∂L∂logθ)
        γ = θ[1]; η = θ[2:1+s]; ζ = θ[2+s:end]
    end
    λ = fft([kernel_shiftinvar(x[i,:],x[1,:],β[k,:],β[l,:],α,γ,η,s)+ζ[k]*(k==l)*(i==1) for i=1:n,k=1:n_β,l=1:n_β],1)
    coeffs_perm = permutedims(hcat(map(i->λ[i,:,:]\y_ft[i,:],1:n)...))
    losses[end] = real.(sum(map(i->logdet(λ[i,:,:]),1:n))+sum(conj.(y_ft).*coeffs_perm))
    γs[end] = γ; ηs[end,:] = η; ζs[end,:] = ζ
    coeffs = real.(ifft(coeffs_perm,1))
end

function mean_post(gp::FastGaussianProcess,x::Vector{Float64},β::Vector{Int64})
    kmat = [kernel_shiftinvar(x,gp._x[i,:],β,gp.β[j,:],gp.α,gp.γ,gp.η,gp.s) for i=1:gp.n,j=1:gp.n_β]
    sum(gp.coeffs.*kmat)
end 
(gp::FastGaussianProcess)(x::Vector{Float64},β::Vector{Int64}) = mean_post(gp,x,β)

function cov_post(gp::FastGaussianProcess,x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64})
    kval = kernel_shiftinvar(x1,x2,β1,β2,gp.α,gp.γ,gp.η,gp.s)
    k1vec = [kernel_shiftinvar(x1,gp._x[i,:],β1,gp.β[j,:],gp.α,gp.γ,gp.η,gp.s) for i=1:gp.n,j=1:gp.n_β]
    k2vec = [kernel_shiftinvar(x2,gp._x[i,:],β2,gp.β[j,:],gp.α,gp.γ,gp.η,gp.s) for i=1:gp.n,j=1:gp.n_β]
    k2vectilde = fft(k2vec,1)
    coeffs = real.(ifft(permutedims(hcat(map(i->gp.λ[i,:,:]\k2vectilde[i,:],1:gp.n)...)),1))
    kval-sum(k1vec.*coeffs)
end

function var_post(gp::FastGaussianProcess,x::Vector{Float64},β::Vector{Int64})
    var = cov_post(gp,x,x,β,β)
    -1e-12<var<0 ? 0. : var
end