module GeneratorsTest

using Test
using ForneyLab
import LinearAlgebra: Diagonal
import ForneyLab: entropiesString, energiesString, freeEnergyString, marginalTableString, inboundString, scheduleString, typeString, vagueString, initializationString, optimizeString, recognitionFactorString, algorithmString, valueString

@testset "typeString" begin
    @test typeString(ForneyLab.SPGaussianMeanPrecisionOutNPP) == "SPGaussianMeanPrecisionOutNPP"
end

@testset "valueString" begin
    @test valueString(1) == "1"
    @test valueString([1]) == "[1]"
    @test valueString(mat(1)) == "mat(1)"
    @test valueString([1 2; 2 1]) == "[1 2; 2 1]"
    @test valueString(Diagonal([1])) == "Diagonal([1])"
end

@testset "inboundString" begin
    f() = 1.0 # Define a function

    # custom inbound
    inbound = Dict(:keyword => false,
                   :g       => f)
    @test inboundString(inbound) == "f"

    inbound = Dict(:g => f)
    @test inboundString(inbound) == "g=f"

    # value inbound
    g = FactorGraph()
    var = Variable()
    inbound = Clamp(var, 1.0)
    inbound.dist_or_msg = Message
    @test inboundString(inbound) == "Message(Univariate, PointMass, m=1.0)"

    # buffer inbound
    g = FactorGraph()
    var = Variable()
    inbound = Clamp(var, 1.0)
    inbound.dist_or_msg = Message
    inbound.buffer_id = :x
    inbound.buffer_index = 1
    @test inboundString(inbound) == "Message(Univariate, PointMass, m=data[:x][1])"

    g = FactorGraph()
    var = Variable()
    inbound = Clamp(var, 1.0)
    inbound.dist_or_msg = Message
    inbound.buffer_id = :x
    inbound.buffer_index = 0
    @test inboundString(inbound) == "Message(Univariate, PointMass, m=data[:x])"

    # marginal
    inbound = MarginalEntry()
    inbound.marginal_id = :x
    @test inboundString(inbound) == "marginals[:x]"
    
    # message
    inbound = ScheduleEntry()
    inbound.schedule_index = 1
    @test inboundString(inbound) == "messages[1]"

    # nothing
    @test inboundString(nothing) == "nothing"
end

@testset "vagueString" begin
    entry = ScheduleEntry()
    entry.family = GaussianMeanPrecision
    entry.dimensionality = ()
    entry_str = vagueString(entry)
    @test entry_str == "vague(GaussianMeanPrecision)"

    entry = ScheduleEntry()
    entry.family = GaussianMeanPrecision
    entry.dimensionality = (1,)
    entry_str = vagueString(entry)
    @test entry_str == "vague(GaussianMeanPrecision, (1,))"
end

@testset "marginalTableString" begin
    inbounds = Vector{ScheduleEntry}(undef, 2)
    inbounds[1] = ScheduleEntry()
    inbounds[1].schedule_index = 1
    inbounds[2] = ScheduleEntry()
    inbounds[2].schedule_index = 2

    marginal_table = Vector{MarginalEntry}(undef, 3)
    marginal_table[1] = MarginalEntry()
    marginal_table[1].marginal_id = :x
    marginal_table[1].marginal_update_rule = Nothing
    marginal_table[1].inbounds = [inbounds[1]]
    marginal_table[2] = MarginalEntry()
    marginal_table[2].marginal_id = :y
    marginal_table[2].marginal_update_rule = ForneyLab.Product
    marginal_table[2].inbounds = inbounds
    marginal_table[3] = MarginalEntry()
    marginal_table[3].marginal_id = :z
    marginal_table[3].marginal_update_rule = ForneyLab.MGaussianMeanPrecisionGGD
    marginal_table[3].inbounds = inbounds

    marginal_table_str = marginalTableString(marginal_table)

    @test occursin("marginals[:x] = messages[1].dist", marginal_table_str)
    @test occursin("marginals[:y] = messages[1].dist * messages[2].dist", marginal_table_str)
    @test occursin("marginals[:z] = ruleMGaussianMeanPrecisionGGD(messages[1], messages[2])", marginal_table_str)
end

@testset "scheduleString" begin
    schedule = Vector{ScheduleEntry}(undef, 2)
    schedule[1] = ScheduleEntry()
    schedule[1].schedule_index = 1
    schedule[1].message_update_rule = Nothing
    schedule[1].inbounds = []
    schedule[2] = ScheduleEntry()
    schedule[2].message_update_rule = ForneyLab.SPGaussianMeanPrecisionOutNPP
    schedule[2].schedule_index = 2
    schedule[2].inbounds = [schedule[1]]

    schedule_str = scheduleString(schedule)
    @test occursin("messages[2] = ruleSPGaussianMeanPrecisionOutNPP(messages[1])", schedule_str)
end

@testset "entropiesString" begin
    inbounds = Vector{MarginalEntry}(undef, 2)
    inbounds[1] = MarginalEntry()
    inbounds[1].marginal_id = :y_x
    inbounds[2] = MarginalEntry()
    inbounds[2].marginal_id = :x

    entropies_vect = [Dict(:conditional => false,
                           :inbounds    => [inbounds[2]]), 
                      Dict(:conditional => true,
                           :inbounds    => inbounds)]

    entropies_str = entropiesString(entropies_vect)

    @test occursin("F -= differentialEntropy(marginals[:x])", entropies_str)
    @test occursin("F -= conditionalDifferentialEntropy(marginals[:y_x], marginals[:x])", entropies_str)
end

@testset "energiesString" begin
    inbound = MarginalEntry()
    inbound.marginal_id = :x
    energies_vect = [Dict(:node     => GaussianMeanPrecision,
                          :inbounds => [inbound])]
    energies_str = energiesString(energies_vect)

    @test occursin("F += averageEnergy(GaussianMeanPrecision, marginals[:x])", energies_str)
end

@testset "initializationString" begin
    algo = Algorithm()
    rf = RecognitionFactor(algo, id=:X)
    rf.initialize = true
    entry = ScheduleEntry()
    entry.schedule_index = 1
    entry.initialize = true
    entry.family = GaussianMeanPrecision
    entry.dimensionality = ()
    rf.schedule = [entry]

    rf_str = initializationString(rf)
    @test occursin("function initX()", rf_str)
    @test occursin("messages[1] = Message(vague(GaussianMeanPrecision))", rf_str)
end

@testset "optimizeString" begin
    algo = Algorithm()
    rf = RecognitionFactor(algo, id=:X)
    rf.optimize = true

    rf_str = optimizeString(rf)
    @test occursin("function optimizeX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=initX()", rf_str)    
end

@testset "recognitionFactorString" begin
    algo = Algorithm()
    rf = RecognitionFactor(algo, id=:X)
    rf.schedule = []
    rf.marginal_table = []

    rf_str = recognitionFactorString(rf)
    @test occursin("function stepX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, 0))", rf_str)
end

@testset "freeEnergyString" begin
    algo = Algorithm()
    free_energy_str = freeEnergyString(algo)
    @test occursin("function freeEnergy(data::Dict, marginals::Dict)", free_energy_str)
end

@testset "algorithmString" begin
    algo = Algorithm()
    algo_str = algorithmString(algo)
    @test occursin("begin", algo_str)
end

end # module