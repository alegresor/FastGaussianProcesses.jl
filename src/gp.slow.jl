mutable struct GaussianProcessRBF
    s::Int64
    n::Int64 
    γ::Float64
    η::Vector{Float64}
    ζ::Vector{Float64}
    β::Matrix{Int64}
    n_β::Int64
    x::Matrix{Float64}
    y::Matrix{Float64}
    ν::Vector{Float64}
    evecs::Matrix{Float64}
    evals::Vector{Float64}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
    optim_steps::Int64
    learnrate::Vector{Float64}
    decayrate::Vector{Float64}
    NEGVARTHRESHOLD::Float64
end

function GaussianProcessRBF(
        x::Matrix{Float64},
        y::Matrix{Float64};
        β::Union{Nothing,Matrix{Int64}} = nothing,
        γ::Float64 = 1.,
        η::Union{Float64,Vector{Float64}} = 1.,
        ζ::Union{Float64,Vector{Float64}} = 1e-12,
        optim_steps::Int64 = 320,
        learnrate_γ::Float64 = 1e-1,
        learnrate_η::Union{Float64,Vector{Float64}} = 1e-1,
        learnrate_ζ::Union{Float64,Vector{Float64}} = 0.,
        decayrate_γ::Float64 = .9,
        decayrate_η::Union{Float64,Vector{Float64}} = .9,
        decayrate_ζ::Union{Float64,Vector{Float64}} = .9,
        verbose::Int64 = 40,
        NEGVARTHRESHOLD::Float64 = -1e-8)
    n,s = size(x)
    if β===nothing β = zeros(Int64,1,s) else @assert size(β,2)==s "β must be a two dimensional matrix of size (n_β,s)" end; @assert β[1,:]==zeros(s); 
    n_β = size(β,1); @assert all(0 .≤ β .≤ 1); @assert all(0 .≤ sum(β,dims=2) .≤ 1)
    @assert size(β,2) == s;   @assert size(y) == (n,n_β) "the number of columns in y must equal the number of rows of β"
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) "η must be a vector of length s or a Float64 which gets copied to each element" end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,n_β) else @assert size(ζ)==(n_β,) "ζ must be a vector of length n_β or a Float64 which gets copied to each element" end
    if typeof(learnrate_η)==Float64 learnrate_η = learnrate_η*ones(Float64,s) else @assert size(learnrate_η)==(s,) "learnrate_η must be a vector of length s or a Float64 which gets copied to each element" end
    if typeof(learnrate_ζ)==Float64 learnrate_ζ = learnrate_ζ*ones(Float64,n_β) else @assert size(learnrate_ζ)==(n_β,) "learnrate_ζ must be a vector of length n_β or a Float64 which gets copied to each element" end
    if typeof(decayrate_η)==Float64 decayrate_η = decayrate_η*ones(Float64,s) else @assert size(decayrate_η)==(s,) "decayrate_η must be a vector of length s or a Float64 which gets copied to each element" end
    if typeof(decayrate_ζ)==Float64 decayrate_ζ = decayrate_ζ*ones(Float64,n_β) else @assert size(decayrate_ζ)==(n_β,) "decayrate_ζ must be a vector of length n_β or a Float64 which gets copied to each element" end
    gp = GaussianProcessRBF(s,n,γ,η,ζ,β,n_β,x,y,Vector{Float64}(undef,n_β*n),Matrix{Float64}(undef,n_β*n,n_β*n),Vector{Float64}(undef,n_β*n),Vector{Float64}(undef,optim_steps+1),Vector{Float64}(undef,optim_steps+1),Matrix{Float64}(undef,optim_steps+1,s),Matrix{Float64}(undef,optim_steps+1,n_β),optim_steps,vcat(learnrate_γ,learnrate_η,learnrate_ζ),vcat(decayrate_γ,decayrate_η,decayrate_ζ),NEGVARTHRESHOLD)
    _train(gp,verbose) 
end
function GaussianProcessRBF(f::Function,seq::IIDU01Seq,n::Int64;kwargs...) 
    x = Next(seq,n)
    y = reshape(vcat([f(x[i,:]) for i=1:n]'...),n,:)
    GaussianProcessRBF(x,y,;kwargs...)
end 
GaussianProcessRBFIIDU01(f::Function,s::Int64,n::Int64;kwargs...) = GaussianProcessRBF(f,IIDU01Seq(s),n;kwargs...)

function _train(gp::GaussianProcessRBF,verbose::Int64)
    verbosebool = verbose > 0
    p = 1+gp.s+gp.n_β    
    k_com = Matrix{Float64}(undef,gp.n_β*gp.n,gp.n_β*gp.n)
    k_nsy = Matrix{Float64}(undef,gp.n_β*gp.n,gp.n_β*gp.n)
    ∂k∂γ = Matrix{Float64}(undef,gp.n_β*gp.n,gp.n_β*gp.n)
    ∂k∂η = Array{Float64}(undef,gp.s,gp.n_β*gp.n,gp.n_β*gp.n)
    ∂k∂ζ = Array{Float64}(undef,gp.n_β,gp.n_β*gp.n,gp.n_β*gp.n)
    ∂k∂θ = Array{Float64}(undef,p,gp.n_β*gp.n,gp.n_β*gp.n)
    ∂L∂θ = Vector{Float64}(undef,p)
    ∂L∂logθ = Vector{Float64}(undef,p)
    θ = Vector{Float64}(undef,p)
    Δ = zeros(Float64,p)
    if verbosebool println("Loss") end
    for step=1:gp.optim_steps+1
        gp.γs[step] = gp.γ; gp.ηs[step,:] .= gp.η; gp.ζs[step,:] .= gp.ζ
        for k=1:gp.n_β,l=1:gp.n_β
            rs,cs = (k-1)*gp.n,(l-1)*gp.n
            for i=1:gp.n,j=1:gp.n k_com[rs+i,cs+j] = rbf_kernel(gp.x[i,:],gp.x[j,:],gp.β[k,:],gp.β[l,:],gp.γ,gp.η) end 
        end
        k_nsy .= k_com; for k=1:gp.n_β k_nsy[(k-1)*gp.n+1:k*gp.n,(k-1)*gp.n+1:k*gp.n] .+= diagm(gp.ζ[k]*ones(Float64,gp.n)) end
        decomp = eigen(Symmetric((k_nsy + k_nsy')./2)); gp.evecs .= decomp.vectors; gp.evals .= max.(decomp.values,1e-32) # k_nsy ≈ evecs*diagm(evals)*evecs'
        gp.ν = gp.evecs*(gp.evecs'*gp.y[:]./gp.evals)
        gp.losses[step] = sum(log.(gp.evals))+gp.y[:]'*gp.ν
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,gp.losses[step]); end
        if step == gp.optim_steps+1 break end 
        ∂k∂γ .= k_com./gp.γ
        for jj=1:gp.s 
            for k=1:gp.n_β,l=1:gp.n_β
                rs,cs = (k-1)*gp.n,(l-1)*gp.n
                for i=1:gp.n,j=1:gp.n ∂k∂η[jj,rs+i,cs+j] = _∂ηl_rbf_kernel(gp.x[i,:],gp.x[j,:],gp.β[k,:],gp.β[l,:],gp.γ,gp.η,jj) end 
            end
        end
        ∂k∂ζ .= 0; for j=1:gp.n_β, k=1:gp.n ∂k∂ζ[j,(j-1)*gp.n+k,(j-1)*gp.n+k] = 1 end
        θ[1] = gp.γ; θ[2:gp.s+1] .= gp.η; θ[gp.s+2:end] .= gp.ζ
        ∂k∂θ[1,:,:] .= ∂k∂γ; 
        ∂k∂θ[2:gp.s+1,:,:] .= ∂k∂η; 
        ∂k∂θ[gp.s+2:end,:,:] .= ∂k∂ζ
        for j=1:p ∂L∂θ[j,:,:] .= sum(map(k->(gp.evecs*(gp.evecs'*∂k∂θ[j,k,:]./gp.evals))[k],1:gp.n_β*gp.n))-gp.ν'*∂k∂θ[j,:,:]*gp.ν end
        ∂L∂logθ .= ∂L∂θ.*θ # since ∂θ∂logθ = θ
        Δ .= gp.decayrate.*Δ.+(1 .- gp.decayrate).*∂L∂logθ.^2
        θ .= exp.(log.(θ)-gp.learnrate.*Δ.^(-1/2).*∂L∂logθ)
        gp.γ = θ[1]; gp.η .= θ[2:1+gp.s]; gp.ζ .= θ[2+gp.s:end]
    end
    gp
end

function rbf_kernel(x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},γ::Float64,η::Vector{Float64})
    @assert all(0 .≤ β1 .≤ 1);@assert all(0 .≤ β2 .≤ 1); @assert all(0 .≤ sum(β1) .≤ 1); @assert all(0 .≤ sum(β2) .≤ 1)
    i1,i2 = findfirst(b->b==1,β1),findfirst(b->b==1,β2)
    k00 = γ*exp(-sum(((x1.-x2).^2)./η))
    if (i1===nothing)&(i2===nothing) # no derivatives 
        return k00
    elseif i2===nothing # ∂x1[i1]
        return -2*k00*(x1[i1]-x2[i1])/η[i1]
    elseif i1===nothing # ∂x2[i2]
        return 2*k00*(x1[i2]-x2[i2])/η[i2]
    elseif i1 ≠ i2 # ∂x1[i1]∂x2[i2] for i1 ≠ i2
        return -4*k00*(x1[i1]-x2[i1])/η[i1]*(x1[i2]-x2[i2])/η[i2]
    elseif β1 == β2 # ∂x1[i]∂x2[i] for i=i1=i2
        i = i1 # = i2
        return -2*k00*(2*(x1[i]-x2[i])^2/η[i]^2 - 1/η[i])
    else @assert false "programming cases error with β1=$β1 and β2=$β2" end 
end 

function _∂ηl_rbf_kernel(x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},γ,η::Vector{Float64},l::Int64)
    @assert all(0 .≤ β1 .≤ 1);@assert all(0 .≤ β2 .≤ 1); @assert all(0 .≤ sum(β1) .≤ 1); @assert all(0 .≤ sum(β2) .≤ 1)
    i1,i2 = findfirst(b->b==1,β1),findfirst(b->b==1,β2)
    k00 = γ*exp(-sum(((x1.-x2).^2)./η))
    if (i1===nothing)&(i2===nothing) # no derivatives 
        return k00*(x1[l]-x2[l])^2/η[l]^2
    elseif i2===nothing # ∂x1[i1]
        if l==i1
            return -2*k00*((x1[l]-x2[l])^3/η[l]^3-(x1[l]-x2[l])/η[l]^2)
        else # l ≠ i1 
            return -2*k00*(x1[i1]-x2[i1])/η[i1]*(x1[l]-x2[l])^2/η[l]^2
        end  
    elseif i1===nothing # ∂x2[i2]
        if l==i2
            return 2*k00*((x1[l]-x2[l])^3/η[l]^3-(x1[l]-x2[l])/η[l]^2)
        else # l ≠ i2 
            return 2*k00*(x1[i2]-x2[i2])/η[i2]*(x1[l]-x2[l])^2/η[l]^2
        end  
    elseif i1 ≠ i2 # ∂x1[i1]∂x2[i2] for i1 ≠ i2
        if l==i1 
            return -4*k00*((x1[l]-x2[l])^3/η[l]^3-(x1[l]-x2[l])/η[l]^2)*(x1[i2]-x2[i2])/η[i2]
        elseif l==i2  
            return -4*k00*((x1[l]-x2[l])^3/η[l]^3-(x1[l]-x2[l])/η[l]^2)*(x1[i1]-x2[i1])/η[i1]
        else # l ∉ {i1,i2}
            return -4*k00*(x1[i1]-x2[i1])/η[i1]*(x1[i2]-x2[i2])/η[i2]*(x1[l]-x2[l])^2/η[l]^2
        end 
    elseif β1 == β2 # ∂x1[i]∂x2[i] for i=i1=i2
        i = i1 # = i2
        if l==i
            return -2*k00*((x1[l]-x2[l])^2/η[l]^2*(2*(x1[l]-x2[l])^2/η[l]^2-1/η[l])+(1/η[l]^2-4*(x1[i]-x2[i])^2/η[i]^3))
        else # l ≠ i
            return -2*k00*(2*(x1[i]-x2[i])^2/η[i]^2-1/η[i])*(x1[l]-x2[l])^2/η[l]^2
        end
    else @assert false "programming cases error with β1=$β1 and β2=$β2" end 
end

function mean_post(gp::GaussianProcessRBF,x::Vector{Float64};β::Union{Nothing,Vector{Int64}}=nothing)
    if β===nothing β = zeros(Int64,gp.s) else @assert size(β)==(gp.s,) "β must be a length s vector" end 
    k1 = [rbf_kernel(x,gp.x[i,:],β,gp.β[k,:],gp.γ,gp.η) for i=1:gp.n,k=1:gp.n_β]
    k1[:]'*gp.ν
end 
(gp::GaussianProcessRBF)(x::Vector{Float64};kwargs...) = mean_post(gp,x;kwargs...)

function cov_post(gp::GaussianProcessRBF,x1::Vector{Float64},x2::Vector{Float64};β1::Union{Nothing,Vector{Int64}}=nothing,β2::Union{Nothing,Vector{Int64}}=nothing)
    if β1===nothing β1 = zeros(Int64,gp.s) else @assert size(β1)==(gp.s,) "β1 must be a length s vector" end 
    if β2===nothing β2 = zeros(Int64,gp.s) else @assert size(β2)==(gp.s,) "β2 must be a length s vector" end 
    kval = rbf_kernel(x1,x2,β1,β2,gp.γ,gp.η)
    k1 = [rbf_kernel(x1,gp.x[i,:],β1,gp.β[k,:],gp.γ,gp.η) for i=1:gp.n,k=1:gp.n_β]
    k2 = [rbf_kernel(x2,gp.x[i,:],β2,gp.β[k,:],gp.γ,gp.η) for i=1:gp.n,k=1:gp.n_β]
    kval - k1[:]'*(gp.evecs*(gp.evecs'k2[:]./gp.evals))
end