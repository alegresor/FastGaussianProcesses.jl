mutable struct GaussianProcess
    NEGVARTHRESHOLD::Float64
end

function GaussianProcess(;NEGVARTHRESHOLD=-1e-12)
    gp = GaussianProcess(NEGVARTHRESHOLD)
    _train(gp,verbose) 
end

function _train(gp::GaussianProcess,verbose::Int64)
    verbosebool = verbose > 0
    p = 1+gp.s+gp.n_β    
    gp
end

function _multidim_kernel(gp::GaussianProcess,x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64})
    k = gp.γ
end 

function mean_post(gp::GaussianProcess,x::Vector{Float64};β::Union{Nothing,Vector{Int64}}=nothing)
    if β===nothing β = zeros(Int64,gp.s) else @assert size(β)==(gp.s,) "β must be a length s vector" end 
end 
(gp::GaussianProcess)(x::Vector{Float64};kwargs...) = mean_post(gp,x;kwargs...)

function cov_post(gp::GaussianProcess,x1::Vector{Float64},x2::Vector{Float64};β1::Union{Nothing,Vector{Int64}}=nothing,β2::Union{Nothing,Vector{Int64}}=nothing)
    if β1===nothing β1 = zeros(Int64,gp.s) else @assert size(β1)==(gp.s,) "β1 must be a length s vector" end 
    if β2===nothing β2 = zeros(Int64,gp.s) else @assert size(β2)==(gp.s,) "β2 must be a length s vector" end 
end