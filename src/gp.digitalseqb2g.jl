function k4sumterm(_x::UInt64,t::Int64)
    total = BigFloat(0.) 
    for a1=0:t-1
        total += (-1)^((_x>>(t-a1-1))&1)/BigFloat(2)^(3*a1)
    end 
    Float64(total)
end 

DSKERNELFUNCSLOW = Dict(
    2 => (β::Int64,x::Float64,_x::UInt64,t::Int64) -> x==0 ? 5/2 : -β*x + 5/2*(1-2. ^(-β)),
    3 => (β::Int64,x::Float64,_x::UInt64,t::Int64) -> x==0 ? 43/18 : β*x^2 - 5*(1-2. ^(-β))*x + 43/18*(1-2. ^(-2*β)),
    4 => (β::Int64,x::Float64,_x::UInt64,t::Int64) -> x==0 ? 701/294 : -2/3*β*x^3 + 5*(1-2. ^(-β))*x^2 - 43/9*(1-2. ^(-2*β))*x + 701/294*(1-2. ^(-3*β)) + β*(1/48*k4sumterm(_x,t)-1/42))

function kernel_digshiftinvar_s1(_x1::UInt64,_x2::UInt64,β1::Int64,β2::Int64,α::Int64,t::Int64)
    _x = _x1 ⊻ _x2
    x = _x * 2. ^(-t)
    β = x == 0 ? 0 : Int64(-floor(log2(x)))
    ktildeidx = α-β1-β2
    (-2)^(β1+β2)*DSKERNELFUNCSLOW[ktildeidx](β,x,_x,t)
end 
kernel_digshiftinvar_s1(_x1::UInt64,x2::Float64,β1::Int64,β2::Int64,α::Int64,t::Int64) = kernel_digshiftinvar_s1(_x1,Float64ToBinary(x2,t),β1,β2,α,t)
kernel_digshiftinvar_s1(x1::Float64,_x2::UInt64,β1::Int64,β2::Int64,α::Int64,t::Int64) = kernel_digshiftinvar_s1(Float64ToBinary(x1,t),_x2,β1,β2,α,t)
kernel_digshiftinvar_s1(x1::Float64,x2::Float64,β1::Int64,β2::Int64,α::Int64,t::Int64) = kernel_digshiftinvar_s1(Float64ToBinary(x1,t),Float64ToBinary(x2,t),β1,β2,α,t)

GaussianProcessDigitalSeqB2G(f::Function,s::Int64,n::Int64;kwargs...) = FastGaussianProcess(f,RandomDigitalShift(DigitalSeqB2G(LinearMatrixScramble(s))),n;kwargs...)
