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
    λ::Union{Array{ComplexF64},Array{Float64}}
    coeffs::Matrix{Float64}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
    ft::Function # ft(a) = V^H*a
    kernel_1s::Function
    optim_steps::Int64
    learningrate::Float64
    decayrate::Float64
end

function FastGaussianProcess(f::Function,seq::Union{LatticeSeqB2,RandomShift,DigitalSeqB2G,RandomDigitalShift},n::Int64;β::Union{Nothing,Matrix{Int64}}=nothing,α::Int64=4,γ::Float64=1.,η::Union{Float64,Vector{Float64}}=1.,ζ::Union{Float64,Vector{Float64}}=1e-16,optim_steps::Int64=100,learningrate::Float64=1e-1,decayrate::Float64=.9,verbose::Int64=10)
    @assert log2(n)%1==0; m = Int64(log2(n))
    if typeof(seq) in [LatticeSeqB2,RandomShift] # lattice sequence
        x = _x = FirstLinear(seq,m)
        @assert 2*α in keys(BERNOULLIPOLYS)
        ft = (a,dims) -> fft(a,dims)/sqrt(n)
        ft_dtype = ComplexF64
        kernel_1s = kernel_shiftinvar_s1
    else # digital sequence
        _x = FirstLinearBinary(seq,m)
        x = BinaryToFloat64(_x,seq)
        ft = (a,dims) -> fwht_natural(a,dims)*sqrt(n)
        ft_dtype = Float64
        kernel_1s = (_x1::Union{UInt64,Float64},_x2::Union{UInt64,Float64},β1::Int64,β2::Int64,α::Int64) -> kernel_digshiftinvar_s1(_x1,_x2,β1,β2,α,seq.t)
    end 
    s = size(x,2)
    if β===nothing β = zeros(Int64,1,s) else @assert size(β,2)==s end
    n_β = size(β,1)
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,n_β) else @assert size(ζ)==(n_β,) end
    y = reshape(vcat([f(x[i,:]) for i=1:n]'...),n,n_β)
    gp = FastGaussianProcess(s,n,α,γ,η,ζ,β,n_β,_x,x,y,Array{ft_dtype}(undef,n,n_β,n_β),Matrix{Float64}(undef,n,n_β),Vector{Float64}(undef,optim_steps+1),Vector{Float64}(undef,optim_steps+1),Matrix{Float64}(undef,optim_steps+1,s),Matrix{Float64}(undef,optim_steps+1,n_β),ft,kernel_1s,optim_steps,learningrate,decayrate)
    _train(gp,verbose) 
end

function _train(gp::FastGaussianProcess,verbose::Int64)
    verbosebool = verbose > 0
    p = 1+gp.s+gp.n_β
    y_ft = conj.(gp.ft(gp.y,1))
    if verbosebool println("QGP Optimization Loss") end
    Δ = zeros(Float64,p)
    kvals_sep = Array{Float64}(undef,gp.n,gp.n_β,gp.n_β,gp.s)
    coeffs_perm = copy(gp.λ)
    for step=1:gp.optim_steps+1
        gp.γs[step] = gp.γ; gp.ηs[step,:] .= gp.η; gp.ζs[step,:] .= gp.ζ
        # make the following more efficient by evaluating K(x_{ij},x_{1j},β_{kj},β_{lj}) for all k,l simultaneously
        kvals_com = gp.γ*ones(Float64,gp.n,gp.n_β,gp.n_β)
        for j=1:gp.s,i=1:gp.n,k=1:gp.n_β,l=1:gp.n_β 
            kvals_sep[i,k,l,j] = gp.kernel_1s(gp._x[i,j],gp._x[1,j],gp.β[k,j],gp.β[l,j],gp.α)
            kvals_com[i,k,l] *= gp.β[k,j]+gp.β[l,j] == 0 ? 1+gp.η[j]*kvals_sep[i,k,l,j] : gp.η[j]*kvals_sep[i,k,l,j]
        end
        kvals_noisy = copy(kvals_com); for k=1:gp.n_β kvals_noisy[1,k,k] += gp.ζ[k] end 
        gp.λ .= sqrt(gp.n).*conj.(gp.ft(kvals_noisy,1)) # block diagonal matrix of eigenvalues, but sparsely stored 3 dim array  
        coeffs_perm = permutedims(hcat(map(i->gp.λ[i,:,:]\y_ft[i,:],1:gp.n)...))
        gp.losses[step] = real.(sum(map(i->logdet(gp.λ[i,:,:]),1:gp.n))+sum(y_ft'*coeffs_perm))
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,gp.losses[step]) end
        if step == gp.optim_steps+1 break end 
        ∂k∂γ = kvals_com./gp.γ
        ∂k∂η = zeros(Float64,gp.s,gp.n,gp.n_β,gp.n_β)
        for j=1:gp.s,i=1:gp.n,k=1:gp.n_β,l=1:gp.n_β
            if kvals_sep[i,k,l,j]==0 ∂k∂η[j,i,k,l] = 0; continue end 
            denom = gp.β[k,j]+gp.β[l,j] == 0 ? 1+gp.η[j]*kvals_sep[i,k,l,j] : gp.η[j]*kvals_sep[i,k,l,j]
            ∂k∂η[j,i,k,l] = kvals_com[i,k,l]*kvals_sep[i,k,l,j]/denom
        end
        ∂k∂ζ = zeros(Float64,gp.n_β,gp.n,gp.n_β,gp.n_β); for k=1:gp.n_β ∂k∂ζ[k,1,k,k] = 1 end
        θ = [gp.γ,gp.η...,gp.ζ...]
        ∂k∂θ = zeros(Float64,p,gp.n,gp.n_β,gp.n_β); ∂k∂θ[1,:,:,:] = ∂k∂γ; ∂k∂θ[2:gp.s+1,:,:,:] = ∂k∂η; ∂k∂θ[gp.s+2:end,:,:,:] = ∂k∂ζ
        ∂ktilde∂θ = sqrt(gp.n).*conj.(gp.ft(∂k∂θ,2))
        ∂L∂θ = real.([sum(map(i->tr(gp.λ[i,:,:]\∂ktilde∂θ[j,i,:,:])-coeffs_perm[i,:]'*∂ktilde∂θ[j,i,:,:]*coeffs_perm[i,:],1:gp.n)) for j=1:p])
        ∂L∂logθ = ∂L∂θ.*θ # since ∂θ∂logθ = θ
        Δ = gp.decayrate*Δ.+(1-gp.decayrate)*∂L∂logθ.^2
        θ = exp.(log.(θ)-gp.learningrate*Δ.^(-1/2).*∂L∂logθ)
        gp.γ = θ[1]; gp.η .= θ[2:1+gp.s]; gp.ζ .= θ[2+gp.s:end]
    end
    gp.coeffs .= real.(gp.ft(coeffs_perm,1))
    gp
end

function mean_post(gp::FastGaussianProcess,x::Union{Vector{Float64},Vector{UInt64}},β::Vector{Int64})
    kvals_sep = [gp.kernel_1s(x[j],gp._x[i,j],β[j],gp.β[k,j],gp.α) for i=1:gp.n,k=1:gp.n_β,j=1:gp.s]
    kmat = gp.γ.*[prod(map(j->β[j]+gp.β[k,j] == 0 ? 1+gp.η[j]*kvals_sep[i,k,j] : gp.η[j]*kvals_sep[i,k,j],1:gp.s))  for i=1:gp.n,k=1:gp.n_β]
    sum(gp.coeffs.*kmat)
end 
(gp::FastGaussianProcess)(x::Union{Vector{Float64},Vector{UInt64}},β::Vector{Int64}) = mean_post(gp,x,β)

function cov_post(gp::FastGaussianProcess,x1::Union{Vector{Float64},Vector{UInt64}},x2::Union{Vector{Float64},Vector{UInt64}},β1::Vector{Int64},β2::Vector{Int64})
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