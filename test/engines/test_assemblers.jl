module AssemblersTest

using Test
using ForneyLab
import ForneyLab: assembleBreaker!, assembleClamp!, assembleAlgorithm!, assembleRecognitionFactor!, assembleSchedule!, assembleInitialization!, assembleMarginalTable!, condense, flatten

@testset "assembleClamp!" begin
    g = FactorGraph()
    nd = Clamp(Variable(), 1.0)
    assembleClamp!(nd, ProbabilityDistribution)    
    @test nd.dist_or_msg == ProbabilityDistribution
end

@testset "assembleBreaker!" begin
    breaker_entry = ScheduleEntry()
    assembleBreaker!(breaker_entry, GaussianMeanPrecision, ())
    @test breaker_entry.family == GaussianMeanPrecision
    @test breaker_entry.initialize == true
    @test breaker_entry.dimensionality == ()

    breaker_entry = ScheduleEntry()
    assembleBreaker!(breaker_entry, Union{Gamma, Wishart}, ())
    @test breaker_entry.family == Gamma
    @test breaker_entry.initialize == true
    @test breaker_entry.dimensionality == ()

    breaker_entry = ScheduleEntry()
    assembleBreaker!(breaker_entry, Union{Gamma, Wishart}, (1,1))
    @test breaker_entry.family == Wishart
    @test breaker_entry.initialize == true
    @test breaker_entry.dimensionality == (1,1)
end

@testset "assembleSchedule!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    algo = Algorithm(rfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(rf)
    @test rf.schedule[3].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test rf.schedule[6].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
end

@testset "assembleInitialization!" begin
    # Expectation propagation
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ Probit(x)
    placeholder(y, :y)
    
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = expectationPropagationSchedule(x)
    algo = Algorithm(rfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(rf)
    assembleInitialization!(rf)
    @test rf.schedule[5].message_update_rule == ForneyLab.EPProbitIn1GP
    @test rf.schedule[3].initialize

    # Nonlinear
    f(z) = z
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ Nonlinear(x, f)
    GaussianMeanPrecision(y, 0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    algo = Algorithm(rfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(rf)
    assembleInitialization!(rf)
    @test rf.schedule[7].message_update_rule == ForneyLab.SPNonlinearIn1GG
    @test rf.schedule[3].initialize

    # Optimize
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    placeholder(x, :x)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    algo = Algorithm(rfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(rf)
    assembleInitialization!(rf)
    @test rf.schedule[3].initialize
end

@testset "assembleMarginalTable!" begin
    # Nothing rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    rf.marginal_table = marginalTable(x)
    algo = Algorithm(rfz)
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleMarginalTable!(rf)
    @test rf.marginal_table[1].marginal_update_rule == Nothing
    @test rf.marginal_table[1].marginal_id == :x
    @test rf.marginal_table[1].inbounds == [rf.schedule[3]]

    # Product rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    rf.marginal_table = marginalTable(x)
    algo = Algorithm(rfz)
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleMarginalTable!(rf)
    @test rf.marginal_table[1].marginal_update_rule == ForneyLab.Product
    @test rf.marginal_table[1].marginal_id == :x
    @test rf.marginal_table[1].inbounds == [rf.schedule[3], rf.schedule[6]] 

    # Marginal rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(x, 1.0)
    GaussianMeanPrecision(y, 0.0, 1.0)
    rfz = RecognitionFactorization([x,y], ids=[:XY])
    rf = rfz.recognition_factors[:XY]
    rf.schedule = variationalSchedule(rf)
    rf.marginal_table = marginalTable(rf)
    algo = Algorithm(rfz)
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleMarginalTable!(rf)
    @test rf.marginal_table[3].marginal_update_rule == ForneyLab.MGaussianMeanPrecisionGGD
    @test rf.marginal_table[3].marginal_id == :y_x
    @test rf.marginal_table[3].inbounds == [rf.schedule[3], rf.schedule[1], g.nodes[:clamp_3]]
end

@testset "assembleRecognitionFactor!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    schedule = sumProductSchedule(x)
    rf.schedule = condense(flatten(schedule))
    rf.marginal_table = marginalTable(x)
    algo = Algorithm()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleRecognitionFactor!(rf)
    @test rf.schedule[1].schedule_index == 1
    @test rf.schedule[1].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test rf.schedule[1].inbounds == [nothing, g.nodes[:clamp_1], g.nodes[:clamp_2]]
    @test rf.schedule[2].schedule_index == 2
    @test rf.schedule[2].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test rf.schedule[2].inbounds == [nothing, g.nodes[:clamp_3], g.nodes[:clamp_4]]
    @test rf.marginal_table[1].marginal_id == :x
    @test rf.marginal_table[1].marginal_update_rule == ForneyLab.Product
    @test rf.marginal_table[1].inbounds == [schedule[3], schedule[6]]
end

@testset "assembleAlgorithm!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    schedule = sumProductSchedule(x)
    rf.schedule = condense(flatten(schedule))
    rf.marginal_table = marginalTable(x)
    algo = Algorithm(rfz)
    assembleAlgorithm!(algo)

    @test length(algo.interface_to_schedule_entry) == 2
    @test length(algo.target_to_marginal_entry) == 1
end

end # module