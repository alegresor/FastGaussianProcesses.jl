BERNOULLIPOLYS = Dict(
    #0 =>  x -> 1,
    1 =>  x -> x-1/2,
    2 =>  x -> x^2-x+1/6,
    3 =>  x -> x^3-3/2*x^2+1/2*x,
    4 =>  x -> x^4-2*x^3+x^2-1/30,
    5 =>  x -> x^5-5/2*x^4+5/3*x^3-1/6*x,
    6 =>  x -> x^6-3*x^5+5/2*x^4-1/2*x^2+1/42,
    7 =>  x -> x^7-7/2*x^6+7/2*x^5-7/6*x^3+1/6*x,
    8 =>  x -> x^8-4*x^7+14/3*x^6-7/3*x^4+2/3*x^2-1/30,
    9 =>  x -> x^9-9/2*x^8+6*x^7-21/5*x^5+2*x^3-3/10*x,
    10 => x -> x^10-5*x^9+15/2*x^8-7*x^6+5*x^4-3/2*x^2+5/66)

function kernel_shiftinvar_s1(x1::Float64,x2::Float64,β1::Int64,β2::Int64,α::Int64)
    bpolyidx = 2*α-β1-β2
    (-1)^(α+1+β2)*prod((bpolyidx+1):(2*α))*BERNOULLIPOLYS[bpolyidx](mod(x1-x2,1))
end 

GaussianProcessLatticeSeqB2(f::Function,s::Int64,n::Int64;kwargs...) = FastGaussianProcess(f,RandomShift(LatticeSeqB2(s)),n;kwargs...)
