module JuliaMessagePassingTest

using Test
using ForneyLab

@testset "Julia messagePassingAlgorithm" begin
    @testset "Read from data" begin
        g = FactorGraph()
        x = Variable()
        m = Variable()
        GaussianMeanVariance(m, constant(0.0), constant(1.0))
        GaussianMeanVariance(x, m, constant(1.0))
        placeholder(x, :x)

        schedule = sumProductSchedule(m)
        marginal_schedule = marginalSchedule(m)

        algo = ForneyLab.messagePassingAlgorithm(schedule, marginal_schedule)

        @test occursin("Array{Message}(undef, 2)", algo)
        @test occursin("messages[1] = ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=0.0), Message(Univariate, PointMass, m=1.0))", algo)
        @test occursin("messages[2] = ruleSPGaussianMeanVarianceMPNP(Message(Univariate, PointMass, m=data[:x]), nothing, Message(Univariate, PointMass, m=1.0))", algo)
        @test occursin("marginals[:variable_2] = messages[1].dist * messages[2].dist", algo)
    end

    @testset "Model with dangling edge" begin
        FactorGraph()
        x = Variable()
        GaussianMeanVariance(x, constant(0.0), constant(1.0))

        schedule = sumProductSchedule(x)
        marginal_schedule = marginalSchedule(x)

        algo = ForneyLab.messagePassingAlgorithm(schedule, marginal_schedule)

        @test occursin("Array{Message}(undef, 1)", algo)
        @test occursin("messages[1] = ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=0.0), Message(Univariate, PointMass, m=1.0))", algo)
        @test occursin("marginals[:variable_1] = messages[1].dist", algo)
    end

    @testset "Estimation of clamped parameter by sum-product rule" begin
        FactorGraph()
        x = Variable()
        GaussianMeanVariance(x, constant([0.0]), constant(mat(1.0)))
        placeholder(x, :x, dims=(1,))

        algo = sumProductAlgorithm(x)

        @test occursin("function optimize!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=init())", algo)
        @test occursin("function init()", algo)
        @test occursin("messages[1] = Message(vague(GaussianMeanVariance, (1,)))", algo)
    end

    @testset "Estimation of clamped parameter by variational rule" begin
        FactorGraph()
        x = Variable()
        GaussianMeanVariance(x, constant([0.0]), constant(mat(1.0)))
        placeholder(x, :x, dims=(1,))

        q = RecognitionFactorization(x, ids=[:X])
        algo = ForneyLab.variationalAlgorithm(q)

        @test occursin("function optimizeX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=initX())", algo)
        @test occursin("function initX()", algo)
        @test occursin("messages[1] = Message(vague(GaussianMeanVariance, (1,)))", algo)
    end
end

@testset "Julia algorithm execution" begin
    # Estimate m from model:
    # m ~ N(0,1)
    # y[i] ~ N(m,1)

    # Resulting algorithm ---
    function step!(marginals::Dict, data::Dict)

    messages = Array{Message}(undef, 6)

    messages[1] = ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=0.0), Message(Univariate, PointMass, m=1.0))
    messages[2] = ruleSPGaussianMeanVarianceMPNP(Message(Univariate, PointMass, m=data[:y][2]), nothing, Message(Univariate, PointMass, m=1.0))
    messages[3] = ruleSPGaussianMeanVarianceMPNP(Message(Univariate, PointMass, m=data[:y][3]), nothing, Message(Univariate, PointMass, m=1.0))
    messages[4] = ruleSPEqualityGaussian(messages[3], messages[2], nothing)
    messages[5] = ruleSPGaussianMeanVarianceMPNP(Message(Univariate, PointMass, m=data[:y][1]), nothing, Message(Univariate, PointMass, m=1.0))
    messages[6] = ruleSPEqualityGaussian(nothing, messages[5], messages[4])

    marginals[:variable_1] = messages[1].dist * messages[6].dist

    end
    # ---

    marginals = Dict()
    data = Dict(:y => [1.0, 2.0, 3.0])
    step!(marginals, data)

    @test marginals[:variable_1] == ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=6.0, w=4.0)
end

end # module