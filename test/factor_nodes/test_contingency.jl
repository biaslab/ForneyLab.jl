module ContingencyTest

using Test
using ForneyLab
import ForneyLab: vague, dims

@testset "Contingency ProbabilityDistribution and Message construction" begin
    @test ProbabilityDistribution(Multivariate, Contingency, p=[0.1 0.4; 0.3 0.2]) == ProbabilityDistribution{Multivariate, Contingency}(Dict(:p=>[0.1 0.4; 0.3 0.2]))
    @test_throws Exception ProbabilityDistribution(Univariate, Contingency)
    @test ProbabilityDistribution(Contingency, p=[0.1 0.4; 0.3 0.2]) == ProbabilityDistribution{Multivariate, Contingency}(Dict(:p=>[0.1 0.4; 0.3 0.2]))
    @test ProbabilityDistribution(Contingency) == ProbabilityDistribution{Multivariate, Contingency}(Dict(:p=>(1/9)*ones(3,3)))
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Contingency, p=[0.1 0.4; 0.3 0.2])) == 2
end

@testset "vague" begin
    @test vague(Contingency) == ProbabilityDistribution(Contingency, p=(1/9)*ones(3,3))
    @test vague(Contingency, (2,3)) == ProbabilityDistribution(Contingency, p=(1/6)*ones(2,3))
end

@testset "differentialEntropy" begin
    @test differentialEntropy(ProbabilityDistribution(Multivariate, Contingency, p=[0.2 0.1 0.7; 0.4 0.3 0.3; 0.1 0.6 0.3])) == 2.78866425274534
end

end # module