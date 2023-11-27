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
kernel_digshiftinvar(x1::Vector{Float64},x2::Vector{Float64},β1::Vector{Int64},β2::Vector{Int64},α::Int64,γ::Float64,η::Vector{Float64},s::Int64) = kernel_digshiftinvar(Float64ToBinary.(x1,53),Float64ToBinary.(x2,53),β1,β2,α,γ,η,s,53)
