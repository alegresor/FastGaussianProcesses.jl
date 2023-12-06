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
    ν::Matrix{Float64}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
    optim_steps::Int64
    learningrate::Float64
    decayrate::Float64
    kernel_1s::Function
    ft::Function # V*a = ft(a),      V^H*a = conj.(ft(conj.(a)))
    NEGVARTHRESHOLD::Float64
end

function FastGaussianProcess(f::Function,seq::Union{LatticeSeqB2,RandomShift,DigitalSeqB2G,RandomDigitalShift},n::Int64;β::Union{Nothing,Matrix{Int64}}=nothing,α::Int64=4,γ::Float64=1.,η::Union{Float64,Vector{Float64}}=1.,ζ::Union{Float64,Vector{Float64}}=1e-16,optim_steps::Int64=320,learningrate::Float64=1e-1,decayrate::Float64=.9,verbose::Int64=40,NEGVARTHRESHOLD::Float64=-1e-8)
    s = typeof(seq) in [LatticeSeqB2,DigitalSeqB2G] ? seq.s : seq.seq.s
    if β===nothing β = zeros(Int64,1,s) else @assert size(β,2)==s "β must be a two dimensional matrix of size (n_β,s)" end
    n_β = size(β,1)
    @assert log2(n)%1==0 "n must be a power of 2."; m = Int64(log2(n))
    if typeof(seq) in [LatticeSeqB2,RandomShift] # lattice sequence
        x = _x = FirstLinear(seq,m)
        @assert 2*α in keys(BERNOULLIPOLYS) "α too large for currently supported kernels."
        @assert 2*α-2*maximum(β) in keys(BERNOULLIPOLYS) "maximum(β) too large for currently supported kernels. Try increasing α if possible."
        ft = (a,dims) -> fft(a,dims)/sqrt(n)
        FTDataType = ComplexF64
        kernel_1s = (x1::Float64,x2::Float64,β1::Int64,β2::Int64) -> kernel_shiftinvar_s1(x1,x2,β1,β2,α)
    else # digital sequence
        _x = FirstLinearBinary(seq,m)
        x = BinaryToFloat64(_x,seq)
        @assert α in keys(DSKERNELFUNCSLOW) "α too large for currently supported kernels"
        @assert α-2*maximum(β) in keys(DSKERNELFUNCSLOW) "maximum(β) too large for currently supported kernels. Try increasing α if possible."
        ft = (a,dims) -> fwht_natural(a,dims)*sqrt(n)
        FTDataType = Float64
        kernel_1s = (_x1::Union{UInt64,Float64},_x2::Union{UInt64,Float64},β1::Int64,β2::Int64) -> kernel_digshiftinvar_s1(_x1,_x2,β1,β2,α,seq.t)
    end 
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) "η must be a vector of length s or a Float64 which gets copied to each element" end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,n_β) else @assert size(ζ)==(n_β,) "ζ must be a vector of length n_β or a Float64 which gets copied to each element" end
    y = reshape(vcat([f(x[i,:]) for i=1:n]'...),n,n_β)
    gp = FastGaussianProcess(s,n,α,γ,η,ζ,β,n_β,_x,x,y,Array{FTDataType}(undef,n,n_β,n_β),Matrix{Float64}(undef,n,n_β),Vector{Float64}(undef,optim_steps+1),Vector{Float64}(undef,optim_steps+1),Matrix{Float64}(undef,optim_steps+1,s),Matrix{Float64}(undef,optim_steps+1,n_β),optim_steps,learningrate,decayrate,kernel_1s,ft,NEGVARTHRESHOLD)
    _train(gp,verbose) 
end

function _train(gp::FastGaussianProcess,verbose::Int64)
    verbosebool = verbose > 0
    p = 1+gp.s+gp.n_β    
    ν_hft = similar(gp.λ,gp.n,gp.n_β)
    k_sep = Array{Float64}(undef,gp.n,gp.n_β,gp.n_β,gp.s)
    k_com = Array{Float64}(undef,gp.n,gp.n_β,gp.n_β)
    k_nsy = Array{Float64}(undef,gp.n,gp.n_β,gp.n_β)
    ∂k∂γ = Array{Float64}(undef,gp.n,gp.n_β,gp.n_β)
    ∂k∂η = Array{Float64}(undef,gp.s,gp.n,gp.n_β,gp.n_β)
    ∂k∂ζ = Array{Float64}(undef,gp.n_β,gp.n,gp.n_β,gp.n_β)
    ∂k∂θ = Array{Float64}(undef,p,gp.n,gp.n_β,gp.n_β)
    ∂λ∂θ = similar(gp.λ,p,gp.n,gp.n_β,gp.n_β)
    ∂L∂θ = Vector{Float64}(undef,p)
    ∂L∂logθ = Vector{Float64}(undef,p)
    θ = Vector{Float64}(undef,p)
    Δ = zeros(Float64,p)
    if verbosebool println("Loss") end
    y_hft = conj.(gp.ft(gp.y,1))
    for step=1:gp.optim_steps+1
        gp.γs[step] = gp.γ; gp.ηs[step,:] .= gp.η; gp.ζs[step,:] .= gp.ζ
        # make the following more efficient by evaluating K(x_{ij},x_{1j},β_{kj},β_{lj}) for all k,l simultaneously
        k_com .= gp.γ
        for i=1:gp.n,k=1:gp.n_β,l=1:gp.n_β,j=1:gp.s
            k_sep[i,k,l,j] = gp.kernel_1s(gp._x[i,j],gp._x[1,j],gp.β[k,j],gp.β[l,j])
            k_com[i,k,l] *= gp.β[k,j]+gp.β[l,j] == 0 ? 1+gp.η[j]*k_sep[i,k,l,j] : gp.η[j]*k_sep[i,k,l,j]
        end
        k_nsy .= k_com; for k=1:gp.n_β k_nsy[1,k,k] += gp.ζ[k] end 
        gp.λ .= sqrt(gp.n).*conj.(gp.ft(k_nsy,1)) # block diagonal matrix of eigenvalues, but sparsely stored 3 dim array 
        for i=1:gp.n ν_hft[i,:] .= gp.λ[i,:,:]\y_hft[i,:] end 
        gp.losses[step] = real.(sum(map(i->logdet(gp.λ[i,:,:])+y_hft[i,:]'*ν_hft[i,:],1:gp.n)))
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,gp.losses[step]) end
        if step == gp.optim_steps+1 break end 
        ∂k∂γ .= k_com./gp.γ
        for j=1:gp.s,i=1:gp.n,k=1:gp.n_β,l=1:gp.n_β
            if k_sep[i,k,l,j]==0 ∂k∂η[j,i,k,l] = 0; continue end 
            denom = gp.β[k,j]+gp.β[l,j] == 0 ? 1+gp.η[j]*k_sep[i,k,l,j] : gp.η[j]*k_sep[i,k,l,j]
            ∂k∂η[j,i,k,l] = k_com[i,k,l]*k_sep[i,k,l,j]/denom
        end
        ∂k∂ζ .= 0; for k=1:gp.n_β ∂k∂ζ[k,1,k,k] = 1 end
        θ[1] = gp.γ; θ[2:gp.s+1] .= gp.η; θ[gp.s+2:end] .= gp.ζ
        ∂k∂θ[1,:,:,:] .= ∂k∂γ; ∂k∂θ[2:gp.s+1,:,:,:] .= ∂k∂η; ∂k∂θ[gp.s+2:end,:,:,:] .= ∂k∂ζ
        ∂λ∂θ .= sqrt(gp.n).*conj.(gp.ft(∂k∂θ,2))
        ∂L∂θ .= real.([sum(map(i->tr(gp.λ[i,:,:]\∂λ∂θ[j,i,:,:])-ν_hft[i,:]'*∂λ∂θ[j,i,:,:]*ν_hft[i,:],1:gp.n)) for j=1:p])
        ∂L∂logθ .= ∂L∂θ.*θ # since ∂θ∂logθ = θ
        Δ .= gp.decayrate*Δ.+(1-gp.decayrate)*∂L∂logθ.^2
        θ .= exp.(log.(θ)-gp.learningrate*Δ.^(-1/2).*∂L∂logθ)
        gp.γ = θ[1]; gp.η .= θ[2:1+gp.s]; gp.ζ .= θ[2+gp.s:end]
    end
    gp.ν .= real.(gp.ft(ν_hft,1))
    gp
end

function _multidim_kernel(gp::FastGaussianProcess,x1::Union{Vector{Float64},Vector{UInt64}},x2::Union{Vector{Float64},Vector{UInt64}},β1::Vector{Int64},β2::Vector{Int64})
    k = gp.γ
    for j=1:gp.s
        kj = gp.kernel_1s(x1[j],x2[j],β1[j],β2[j])
        k *= β1[j]+β2[j] == 0 ? 1+gp.η[j]*kj : gp.η[j]*kj
    end 
    k
end 

function mean_post(gp::FastGaussianProcess,x::Union{Vector{Float64},Vector{UInt64}};β::Union{Nothing,Vector{Int64}}=nothing)
    if β===nothing β = zeros(Int64,gp.s) else @assert size(β)==(gp.s,) "β must be a length s vector" end 
    kmat = [_multidim_kernel(gp,x,gp._x[i,:],β,gp.β[k,:]) for i=1:gp.n,k=1:gp.n_β] 
    sum(gp.ν.*kmat)
end 
(gp::FastGaussianProcess)(x::Union{Vector{Float64},Vector{UInt64}};kwargs...) = mean_post(gp,x;kwargs...)

function cov_post(gp::FastGaussianProcess,x1::Union{Vector{Float64},Vector{UInt64}},x2::Union{Vector{Float64},Vector{UInt64}};β1::Union{Nothing,Vector{Int64}}=nothing,β2::Union{Nothing,Vector{Int64}}=nothing)
    if β1===nothing β1 = zeros(Int64,gp.s) else @assert size(β1)==(gp.s,) "β1 must be a length s vector" end 
    if β2===nothing β2 = zeros(Int64,gp.s) else @assert size(β2)==(gp.s,) "β2 must be a length s vector" end 
    kval = _multidim_kernel(gp,x1,x2,β1,β2)
    k1 = [_multidim_kernel(gp,x1,gp._x[i,:],β1,gp.β[k,:]) for i=1:gp.n,k=1:gp.n_β] 
    k2 = [_multidim_kernel(gp,x2,gp._x[i,:],β2,gp.β[k,:]) for i=1:gp.n,k=1:gp.n_β] 
    k1_hft = conj.(gp.ft(k1,1)); k2_hft = conj.(gp.ft(k2,1))
    ν = similar(gp.λ,gp.n,gp.n_β); for i=1:gp.n ν[i,:] .= gp.λ[i,:,:]\k2_hft[i,:] end
    kval-real.(sum(conj.(k1_hft).*ν))
end

function var_post(gp::FastGaussianProcess,x::Vector{Float64};β::Union{Nothing,Vector{Int64}}=nothing)
    var = cov_post(gp,x,x;β1=β,β2=β)
    @assert var > gp.NEGVARTHRESHOLD "variance less than NEGVARTHRESHOLD = $NEGVARTHRESHOLD"
    gp.NEGVARTHRESHOLD<var<0 ? 0. : var
end