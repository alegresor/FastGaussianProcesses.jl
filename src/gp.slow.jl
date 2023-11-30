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
    L::LowerTriangular{Float64, Matrix{Float64}}
    losses::Vector{Float64}
    γs::Vector{Float64}
    ηs::Matrix{Float64}
    ζs::Matrix{Float64}
    optim_steps::Int64
    learningrate::Float64
    decayrate::Float64
    NEGVARTHRESHOLD::Float64
end

function GaussianProcessRBF(x::Matrix{Float64},y::Matrix{Float64};β::Union{Nothing,Matrix{Int64}}=nothing,γ::Float64=1.,η::Union{Float64,Vector{Float64}}=1.,ζ::Union{Float64,Vector{Float64}}=1e-16,optim_steps::Int64=320,learningrate::Float64=1e-1,decayrate::Float64=.9,verbose::Int64=40,NEGVARTHRESHOLD::Float64=-1e-12)
    n,s = size(x)
    if β===nothing β = zeros(Int64,1,s) else @assert size(β,2)==s "β must be a two dimensional matrix of size (n_β,s)" end; @assert β[1,:]==zeros(s); 
    n_β = size(β,1); @assert all(0 .≤ β .≤ 1); @assert all(0 .≤ sum(β,dims=2) .≤ 1)
    @assert size(β,2) == s;   @assert size(y) == (n,n_β) "the number of columns in y must equal the number of rows of β"
    if typeof(η)==Float64 η = η*ones(Float64,s) else @assert size(η)==(s,) "η must be a vector of length s or a Float64 which gets copied to each element" end 
    if typeof(ζ)==Float64 ζ = ζ*ones(Float64,n_β) else @assert size(ζ)==(n_β,) "ζ must be a vector of length n_β or a Float64 which gets copied to each element" end
    gp = GaussianProcessRBF(s,n,γ,η,ζ,β,n_β,x,y,Vector{Float64}(undef,n_β*n),LowerTriangular(Matrix{Float64}(undef,n_β*n,n_β*n)),Vector{Float64}(undef,optim_steps+1),Vector{Float64}(undef,optim_steps+1),Matrix{Float64}(undef,optim_steps+1,s),Matrix{Float64}(undef,optim_steps+1,n_β),optim_steps,learningrate,decayrate,NEGVARTHRESHOLD)
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
    θ = Vector{Float64}(undef,p)
    Δ = zeros(Float64,p)
    if verbosebool println("Loss") end
    for step=1:gp.optim_steps+1
        gp.γs[step] = gp.γ; gp.ηs[step,:] .= gp.η; gp.ζs[step,:] .= gp.ζ
        for k=1:gp.n_β,l=1:gp.n_β
            rs,cs = (k-1)*gp.n,(l-1)*gp.n
            for i=1:gp.n,j=1:gp.n k_com[rs+i,cs+j] = _rbf_kernel(gp.x[i,:],gp.x[j,:],gp.β[k,:],gp.β[l,:],gp.γ,gp.η) end 
        end
        k_nsy .= k_com; for k=1:gp.n_β k_nsy[(k-1)*gp.n+1:k*gp.n,(k-1)*gp.n+1:k*gp.n] .+= diagm(gp.ζ[k]*ones(Float64,gp.n)) end
        gp.L .= cholesky(k_nsy,NoPivot()).L # k_nsy = gp.L*gp.L'
        gp.ν = gp.L'\(gp.L\gp.y[:]) # k_nsy^{-1}y  solves  L(L'ν) = y
        gp.losses[step] = 2*logdet(gp.L)+gp.y[:]'*gp.ν
        if verbosebool && step%verbose==0 @printf("\tstep %-7i %.1e\n",step,gp.losses[step]) end
        if step == gp.optim_steps+1 break end 
        @assert false "hyperparameter optimization not implemented"
    end
    gp
end

function _rbf_kernel(x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},γ::Float64,η::Vector{Float64})
    @assert all(0 .≤ β1 .≤ 1);@assert all(0 .≤ β2 .≤ 1); @assert all(0 .≤ sum(β1) .≤ 1); @assert all(0 .≤ sum(β2) .≤ 1);
    i1,i2 = findfirst(b->b==1,β1),findfirst(b->b==1,β2)
    if (i1===nothing)&(i2===nothing) # no derivatives 
        return γ*exp(-sum(((x1.-x2).^2)./η))
    elseif i2===nothing # ∂x[i1]
        return γ*exp(-sum(((x1.-x2).^2)./η))*(-2*(x1[i1]-x2[i1])/η[i1])
    elseif i1===nothing # ∂x[i2]
        return γ*exp(-sum(((x1.-x2).^2)./η))*(2*(x1[i2]-x2[i2])/η[i2])
    elseif i1 ≠ i2 # ∂x1[i1]∂x2[i2] for i1 ≠ i2
        return γ*exp(-sum(((x1.-x2).^2)./η))*(-2*(x1[i1]-x2[i1])/η[i1])*(2*(x1[i2]-x2[i2])/η[i2])
    elseif β1 == β2 # ∂x1[i]∂x2[i] for i=i1=i2
        i = i1 # =i2
        γ*exp(-sum(((x1.-x2).^2)./η))*(-4*((x1[i]-x2[i])/η[i])^2 + 2/η[i])
    else @assert false "programming cases error with β1=$β1 and β2=$β2" end 
end 

function mean_post(gp::GaussianProcessRBF,x::Vector{Float64};β::Union{Nothing,Vector{Int64}}=nothing)
    if β===nothing β = zeros(Int64,gp.s) else @assert size(β)==(gp.s,) "β must be a length s vector" end 
    k1 = [_rbf_kernel(x,gp.x[i,:],β,gp.β[k,:],gp.γ,gp.η) for i=1:gp.n,k=1:gp.n_β]
    k1[:]'*gp.ν
end 
(gp::GaussianProcessRBF)(x::Vector{Float64};kwargs...) = mean_post(gp,x;kwargs...)

function cov_post(gp::GaussianProcessRBF,x1::Vector{Float64},x2::Vector{Float64};β1::Union{Nothing,Vector{Int64}}=nothing,β2::Union{Nothing,Vector{Int64}}=nothing)
    if β1===nothing β1 = zeros(Int64,gp.s) else @assert size(β1)==(gp.s,) "β1 must be a length s vector" end 
    if β2===nothing β2 = zeros(Int64,gp.s) else @assert size(β2)==(gp.s,) "β2 must be a length s vector" end 
    kval = _rbf_kernel(x1,x2,β1,β2,gp.γ,gp.η)
    k1 = [_rbf_kernel(x1,gp.x[i,:],β1,gp.β[k,:],gp.γ,gp.η) for i=1:gp.n,k=1:gp.n_β]
    k2 = [_rbf_kernel(x2,gp.x[i,:],β2,gp.β[k,:],gp.γ,gp.η) for i=1:gp.n,k=1:gp.n_β]
    kval - k1[:]'*(gp.L'\(gp.L\k2[:]))
end