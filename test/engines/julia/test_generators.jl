module GeneratorsTest

using Test
using ForneyLab
using LinearAlgebra: Diagonal
using ForneyLab: entropiesSourceCode, energiesSourceCode, freeEnergySourceCode, marginalTableSourceCode, inboundSourceCode, scheduleSourceCode, removePrefix, vagueSourceCode, initializationSourceCode, optimizeSourceCode, algorithmSourceCode, valueSourceCode, posteriorFactorSourceCode

@testset "removePrefix" begin
    @test removePrefix(ForneyLab.SPGaussianMeanPrecisionOutNPP) == "SPGaussianMeanPrecisionOutNPP"
end

@testset "valueSourceCode" begin
    @test valueSourceCode(1) == "1"
    @test valueSourceCode([1]) == "[1]"
    @test valueSourceCode(mat(1)) == "mat(1)"
    @test valueSourceCode(reshape([1,2,3], 3, 1)) == "reshape([1, 2, 3], (3, 1))"
    @test valueSourceCode([1 2; 2 1]) == "[1 2; 2 1]"
    @test valueSourceCode(Diagonal([1])) == "Diagonal([1])"
end

f() = 1.0 # Define a function

@testset "inboundSourceCode" begin
    # custom inbound
    inbound = Dict(:keyword => false,
                   :inx     => 1)
    @test inboundSourceCode(inbound) == "1"

    inbound = Dict(:keyword => false,
                   :g       => f)
    @test inboundSourceCode(inbound) == "f"

    inbound = Dict(:g => f)
    @test inboundSourceCode(inbound) == "g=f"

    # value inbound
    g = FactorGraph()
    var = Variable()
    inbound = Clamp(var, 1.0)
    inbound.dist_or_msg = Message
    @test inboundSourceCode(inbound) == "Message(Univariate, PointMass, m=1.0)"

    # buffer inbound
    g = FactorGraph()
    var = Variable()
    inbound = Clamp(var, 1.0)
    inbound.dist_or_msg = Message
    inbound.buffer_id = :x
    inbound.buffer_index = 1
    @test inboundSourceCode(inbound) == "Message(Univariate, PointMass, m=data[:x][1])"

    g = FactorGraph()
    var = Variable()
    inbound = Clamp(var, 1.0)
    inbound.dist_or_msg = Message
    inbound.buffer_id = :x
    inbound.buffer_index = 0
    @test inboundSourceCode(inbound) == "Message(Univariate, PointMass, m=data[:x])"

    # marginal
    inbound = MarginalEntry()
    inbound.marginal_id = :x
    @test inboundSourceCode(inbound) == "marginals[:x]"
    
    # message
    inbound = ScheduleEntry()
    inbound.schedule_index = 1
    @test inboundSourceCode(inbound) == "messages[1]"

    # nothing
    @test inboundSourceCode(nothing) == "nothing"
end

@testset "vagueSourceCode" begin
    entry = ScheduleEntry()
    entry.family = GaussianMeanPrecision
    entry.dimensionality = ()
    entry_code = vagueSourceCode(entry)
    @test entry_code == "vague(GaussianMeanPrecision)"

    entry = ScheduleEntry()
    entry.family = GaussianMeanPrecision
    entry.dimensionality = (1,)
    entry_code = vagueSourceCode(entry)
    @test entry_code == "vague(GaussianMeanPrecision, (1,))"
end

@testset "marginalTableSourceCode" begin
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

    marginal_table_code = marginalTableSourceCode(marginal_table)

    @test occursin("marginals[:x] = messages[1].dist", marginal_table_code)
    @test occursin("marginals[:y] = messages[1].dist * messages[2].dist", marginal_table_code)
    @test occursin("marginals[:z] = ruleMGaussianMeanPrecisionGGD(messages[1], messages[2])", marginal_table_code)
end

@testset "scheduleSourceCode" begin
    schedule = Vector{ScheduleEntry}(undef, 2)
    schedule[1] = ScheduleEntry()
    schedule[1].schedule_index = 1
    schedule[1].message_update_rule = Nothing
    schedule[1].inbounds = []
    schedule[2] = ScheduleEntry()
    schedule[2].message_update_rule = ForneyLab.SPGaussianMeanPrecisionOutNPP
    schedule[2].schedule_index = 2
    schedule[2].inbounds = [schedule[1]]

    schedule_code = scheduleSourceCode(schedule)
    @test occursin("messages[2] = ruleSPGaussianMeanPrecisionOutNPP(messages[1])", schedule_code)
end

@testset "entropiesSourceCode" begin
    inbounds = Vector{MarginalEntry}(undef, 2)
    inbounds[1] = MarginalEntry()
    inbounds[1].marginal_id = :y_x
    inbounds[2] = MarginalEntry()
    inbounds[2].marginal_id = :x

    entropies_vect = [Dict(:counting_number => -1,
                           :inbound         => inbounds[1]), 
                      Dict(:counting_number => 2,
                           :inbound         => inbounds[2])]

    entropies_code = entropiesSourceCode(entropies_vect)

    @test occursin("F -= -1*differentialEntropy(marginals[:y_x])", entropies_code)
    @test occursin("F -= 2*differentialEntropy(marginals[:x])", entropies_code)
end

@testset "energiesSourceCode" begin
    inbound = MarginalEntry()
    inbound.marginal_id = :x
    energies_vect = [Dict(:counting_number => 1,
                          :node            => GaussianMeanPrecision,
                          :inbounds        => [inbound])]
    energies_code = energiesSourceCode(energies_vect)

    @test occursin("F += averageEnergy(GaussianMeanPrecision, marginals[:x])", energies_code)
end

@testset "initializationSourceCode" begin
    g = FactorGraph()
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(pfz, id=:X)
    algo = InferenceAlgorithm(pfz)
    pf.algorithm_id = algo.id
    pf.initialize = true
    entry = ScheduleEntry()
    entry.schedule_index = 1
    entry.initialize = true
    entry.family = GaussianMeanPrecision
    entry.dimensionality = ()
    pf.schedule = [entry]

    pf_code = initializationSourceCode(pf)
    @test occursin("function initX()", pf_code)
    @test occursin("messages[1] = Message(vague(GaussianMeanPrecision))", pf_code)
end

@testset "optimizeSourceCode" begin
    g = FactorGraph()
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(pfz, id=:X)
    algo = InferenceAlgorithm(pfz)
    pf.algorithm_id = algo.id
    pf.optimize = true

    pf_code = optimizeSourceCode(pf)
    @test occursin("function optimizeX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=initX()", pf_code)    
end

@testset "posteriorFactorSourceCode" begin
    g = FactorGraph()
    pfz = PosteriorFactorization() 
    pf = PosteriorFactor(pfz, id=:X)
    algo = InferenceAlgorithm(pfz)
    pf.algorithm_id = algo.id
    pf.schedule = []
    pf.marginal_table = []

    pf_code = posteriorFactorSourceCode(pf)
    @test occursin("function stepX!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, 0))", pf_code)
end

@testset "algorithmSourceCode" begin
    algo = InferenceAlgorithm()
    algo_code = algorithmSourceCode(algo)
    @test occursin("begin", algo_code)
end

@testset "freeEnergySourceCode" begin
    algo = InferenceAlgorithm()
    algo.posterior_factorization.free_energy_flag = true
    free_energy_code = freeEnergySourceCode(algo)
    @test occursin("function freeEnergy(data::Dict, marginals::Dict)", free_energy_code)
end

end # module