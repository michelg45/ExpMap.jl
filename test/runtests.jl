
using  LinearAlgebra
using ExpMap
using Test

@testset "ExpMap.jl" begin
    include("tests_ExpMap.jl")
    res = tests_ExpMap()
    @test res < 1e-12
end
