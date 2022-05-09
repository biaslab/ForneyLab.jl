module ContingencyTest

using Test
using ForneyLab
using ForneyLab: vague, dims

@testset "Contingency Distribution and Message construction" begin
    @test Distribution(Multivariate, Contingency, p=[0.1 0.4; 0.3 0.2]) == Distribution{Multivariate, Contingency}(Dict(:p=>[0.1 0.4; 0.3 0.2]))
    @test_throws Exception Distribution(Univariate, Contingency)
    @test Distribution(Contingency, p=[0.1 0.4; 0.3 0.2]) == Distribution{Multivariate, Contingency}(Dict(:p=>[0.1 0.4; 0.3 0.2]))
    @test Distribution(Contingency) == Distribution{Multivariate, Contingency}(Dict(:p=>(1/9)*ones(3,3)))
end

@testset "dims" begin
    @test dims(Distribution(Contingency, p=[0.1 0.4; 0.3 0.2])) == (2,)
end

@testset "vague" begin
    @test vague(Contingency) == Distribution(Contingency, p=(1/9)*ones(3,3))
    @test vague(Contingency, (2,3)) == Distribution(Contingency, p=(1/6)*ones(2,3))
end

@testset "differentialEntropy" begin
    @test differentialEntropy(Distribution(Multivariate, Contingency, p=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == 2.78866425274534
    @test differentialEntropy(Distribution(Multivariate, Contingency, p=[0.3*[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3], 0.7*[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3]])) == 4.621257158910021
end

end # module