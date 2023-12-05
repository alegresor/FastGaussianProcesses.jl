using Test
using FastGaussianProcesses

@testset "Util" begin
    @test RationalToBinary(1//2,4) == 8
    @test RationalToBinary(1//2,5) == 16
    @test RationalToBinary(5//8,5) == 20
    @test Float64ToBinary(1/2,4) == 8
    @test Float64ToBinary(1/2,5) == 16
    @test Float64ToBinary(5/8,5) == 20
    end

