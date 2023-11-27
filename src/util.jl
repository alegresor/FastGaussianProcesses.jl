function RationalToBinary(r::Rational{Int64},t::Int64)
    b = UInt64(0)
    for k=1:t
        frac = 1//UInt64(2)^k
        if r >= frac
            r -= frac
            b += UInt64(2)^(t-k)
        end
    end
    return b
end

Float64ToBinary(x::Float64,t::Int64) = RationalToBinary(convert(Rational{Int64},x),t)