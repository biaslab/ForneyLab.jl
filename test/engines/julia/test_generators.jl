module GeneratorsTest

using Test
using ForneyLab
import ForneyLab: entropiesString, energiesString, freeEnergyString, marginalScheduleString, inboundString, scheduleString, typeString, vagueString, initializationString, optimizeString, recognitionFactorString, algorithmString

@testset "typeString" begin
    @test typeString(ForneyLab.SPGaussianMeanPrecisionOutNPP) == "SPGaussianMeanPrecisionOutNPP"
end

@testset "inboundString" begin
    f() = 1.0 # Define a function

    # customString
    inbound_dict = Dict(:keyword => false,
                        :g       => f)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "f"

    inbound_dict = Dict(:g => f)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "g=f"

    # valueString
    inbound_dict = Dict(:dist_or_msg  => Message,
                        :variate_type => Univariate,
                        :value        => 1.0)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "Message(Univariate, PointMass, m=1.0)"

    inbound_dict = Dict(:dist_or_msg  => ProbabilityDistribution,
                        :variate_type => MatrixVariate,
                        :value        => mat(1.0))
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "ProbabilityDistribution(MatrixVariate, PointMass, m=mat(1.0))"

    # bufferString
    inbound_dict = Dict(:dist_or_msg  => Message,
                        :variate_type => Univariate,
                        :buffer_id    => :x)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "Message(Univariate, PointMass, m=data[:x])"

    inbound_dict = Dict(:dist_or_msg  => ProbabilityDistribution,
                        :variate_type => MatrixVariate,
                        :buffer_id    => :x,
                        :buffer_index => 1)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "ProbabilityDistribution(MatrixVariate, PointMass, m=data[:x][1])"

    # marginal
    inbound_dict = Dict(:marginal_id => :x)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "marginals[:x]"
    
    # message
    inbound_dict = Dict(:schedule_index => 1)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "messages[1]"

    # nothing
    inbound_dict = Dict(:nothing => true)
    inbound_str = inboundString(inbound_dict)
    @test inbound_str == "nothing"
end

@testset "vagueString" begin
    entry_dict = Dict(:family         => GaussianMeanPrecision,
                        :dimensionality => ())
    entry_str = vagueString(entry_dict)
    @test entry_str == "vague(GaussianMeanPrecision)"

    entry_dict = Dict(:family         => GaussianMeanPrecision,
                        :dimensionality => (1,))
    entry_str = vagueString(entry_dict)
    @test entry_str == "vague(GaussianMeanPrecision, (1,))"
end

@testset "marginalScheduleString" begin
    marginal_schedule = [Dict(:marginal_id          => :x,
                              :marginal_update_rule => Nothing,
                              :inbounds             => [Dict(:schedule_index => 1)]),
                         Dict(:marginal_id          => :y,
                              :marginal_update_rule => ForneyLab.Product,
                              :inbounds             => [Dict(:schedule_index => 1),
                                                        Dict(:schedule_index => 2)]),
                         Dict(:marginal_id          => :z,
                              :marginal_update_rule => ForneyLab.MGaussianMeanPrecisionGGD,
                              :inbounds             => [Dict(:schedule_index => 1),
                                                        Dict(:schedule_index => 2)])]
    marginal_schedule_str = marginalScheduleString(marginal_schedule)

    @test occursin("marginals[:x] = messages[1].dist", marginal_schedule_str)
    @test occursin("marginals[:y] = messages[1].dist * messages[2].dist", marginal_schedule_str)
    @test occursin("marginals[:z] = ruleMGaussianMeanPrecisionGGD(messages[1], messages[2])", marginal_schedule_str)
end

@testset "scheduleString" begin
    schedule = [Dict(:message_update_rule => ForneyLab.SPGaussianMeanPrecisionOutNPP,
                     :schedule_index      => 2,
                     :inbounds            => [Dict(:schedule_index => 1)]),
                Dict(:message_update_rule => ForneyLab.SPGaussianMeanPrecisionOutNPP,
                     :schedule_index      => 3,
                     :inbounds            => [Dict(:schedule_index => 1)])]
    schedule_str = scheduleString(schedule)
    @test occursin("messages[2] = ruleSPGaussianMeanPrecisionOutNPP(messages[1])", schedule_str)
    @test occursin("messages[3] = ruleSPGaussianMeanPrecisionOutNPP(messages[1])", schedule_str)
end

@testset "entropiesString" begin
    entropies_vect = [Dict(:conditional => false,
                           :inbounds    => [Dict(:marginal_id => :x)]), 
                      Dict(:conditional => true,
                           :inbounds    => [Dict(:marginal_id => :y_x),
                                            Dict(:marginal_id => :x)])]
    entropies_str = entropiesString(entropies_vect)

    @test occursin("F -= differentialEntropy(marginals[:x])", entropies_str)
    @test occursin("F -= conditionalDifferentialEntropy(marginals[:y_x], marginals[:x])", entropies_str)
end

@testset "energiesString" begin
    energies_vect = [Dict(:node     => GaussianMeanPrecision,
                          :inbounds => [Dict(:marginal_id => :x)])]
    energies_str = energiesString(energies_vect)

    @test occursin("F += averageEnergy(GaussianMeanPrecision, marginals[:x])", energies_str)
end

@testset "initializationString" begin
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz, id=:X)
    rf.initialize = true
    rf.schedule = [Dict(:schedule_index => 1,
                        :initialize     => true,
                        :family         => GaussianMeanPrecision,
                        :dimensionality => ())]
    rf_str = initializationString(rf)
    @test occursin("function initX()", rf_str)
    @test occursin("messages[1] = Message(vague(GaussianMeanPrecision))", rf_str)
end

@testset "optimizeString" begin
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz, id=:X)
    rf.optimize = true
    rf_str = optimizeString(rf)
    @test occursin("function optimizeX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=initX()", rf_str)    
end

@testset "recognitionFactorString" begin
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz, id=:X)
    rf.schedule = []
    rf.marginal_schedule = []
    rf_str = recognitionFactorString(rf)
    @test occursin("function stepX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, 0))", rf_str)
end

@testset "freeEnergyString" begin
    rfz = RecognitionFactorization()
    free_energy_str = freeEnergyString(rfz)

    @test occursin("function freeEnergy(data::Dict, marginals::Dict)", free_energy_str)
end

@testset "algorithmString" begin
    rfz = RecognitionFactorization()
    algo_str = algorithmString(rfz)
    @test occursin("begin", algo_str)
end

end # module