module AssemblersTest

using Test
using ForneyLab
import ForneyLab: assembleBreaker!, assembleClamp!, assembleAlgorithm!, assembleSchedule!, assembleInitialization!, assembleMarginalTable!, condense, flatten

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
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    rf.schedule = sumProductSchedule(x)
    rf.target_to_marginal_entry = Dict()
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
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
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    rf.schedule = expectationPropagationSchedule(x)
    rf.target_to_marginal_entry = Dict()
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
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
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    rf.schedule = sumProductSchedule(x)
    rf.target_to_marginal_entry = Dict()
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    assembleSchedule!(rf)
    assembleInitialization!(rf)
    @test rf.schedule[7].message_update_rule == ForneyLab.SPNonlinearIn1GG
    @test rf.schedule[3].initialize

    # Optimize
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    placeholder(x, :x)
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    rf.schedule = sumProductSchedule(x)
    rf.target_to_marginal_entry = Dict()
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    assembleSchedule!(rf)
    assembleInitialization!(rf)
    @test rf.schedule[3].initialize
end

@testset "assembleMarginalTable!" begin
    # Nothing rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    rf.schedule = sumProductSchedule(x)
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    rf.marginal_table = marginalTable(x)
    rf.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(rf.marginal_table)
    assembleMarginalTable!(rf)
    @test rf.marginal_table[1].marginal_update_rule == Nothing
    @test rf.marginal_table[1].marginal_id == :x
    @test rf.marginal_table[1].inbounds == [rf.schedule[3]]

    # Product rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    rf.schedule = sumProductSchedule(x)
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    rf.marginal_table = marginalTable(x)
    rf.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(rf.marginal_table)
    assembleMarginalTable!(rf)
    @test rf.marginal_table[1].marginal_update_rule == ForneyLab.Product
    @test rf.marginal_table[1].marginal_id == :x
    @test rf.marginal_table[1].inbounds == [rf.schedule[3], rf.schedule[6]] 

    # Marginal rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(x, 1.0)
    GaussianMeanPrecision(y, 0.0, 1.0)
    algo = Algorithm([x,y], ids=[:XY])
    rf = algo.recognition_factors[:XY]
    rf.schedule = variationalSchedule(rf)
    rf.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    rf.marginal_table = marginalTable(rf)
    rf.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(rf.marginal_table)
    assembleMarginalTable!(rf)
    @test rf.marginal_table[3].marginal_update_rule == ForneyLab.MGaussianMeanPrecisionGGD
    @test rf.marginal_table[3].marginal_id == :y_x
    @test rf.marginal_table[3].inbounds == [rf.schedule[3], rf.schedule[1], g.nodes[:clamp_3]]
end

@testset "assembleAlgorithm!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    algo = Algorithm()
    rf = RecognitionFactor(algo)
    schedule = sumProductSchedule(x)
    rf.schedule = condense(flatten(schedule))
    rf.marginal_table = marginalTable(x)
    assembleAlgorithm!(rf)
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

end # module