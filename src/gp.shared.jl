function var_post(gp::Union{FastGaussianProcess,GaussianProcessRBF},x::Union{Vector{Float64},Vector{UInt64}};β::Union{Nothing,Vector{Int64}}=nothing)
    var = cov_post(gp,x,x;β1=β,β2=β)
    @assert var > gp.NEGVARTHRESHOLD "variance less than NEGVARTHRESHOLD = $NEGVARTHRESHOLD"
    gp.NEGVARTHRESHOLD<var<0 ? 0. : var
end