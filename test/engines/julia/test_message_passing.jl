module JuliaMessagePassingTest

using Test
using ForneyLab
import ForneyLab: assembleAlgorithm!, assembleSchedule!, assembleInitialization!, assembleBreaker!, assembleMarginalSchedule!, assembleClamp!, condense, flatten

@testset "assembleClamp!" begin
    g = FactorGraph()
    nd = Clamp(Variable(), 1.0)
    assembleClamp!(nd, ProbabilityDistribution)    
    @test nd.dist_or_msg == ProbabilityDistribution
end

@testset "assembleMarginalSchedule!" begin
    # Nothing rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(schedule)
    marginal_schedule = marginalSchedule(x)
    target_to_marginal_entry = ForneyLab.targetToMarginalEntry(marginal_schedule)
    assembleMarginalSchedule!(marginal_schedule, interface_to_schedule_entry, target_to_marginal_entry)
    @test marginal_schedule[1].marginal_update_rule == Nothing
    @test marginal_schedule[1].marginal_id == :x
    @test marginal_schedule[1].inbounds == [schedule[3]]

    # Product rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(schedule)
    marginal_schedule = marginalSchedule(x)
    target_to_marginal_entry = ForneyLab.targetToMarginalEntry(marginal_schedule)
    assembleMarginalSchedule!(marginal_schedule, interface_to_schedule_entry, target_to_marginal_entry)
    @test marginal_schedule[1].marginal_update_rule == ForneyLab.Product
    @test marginal_schedule[1].marginal_id == :x
    @test marginal_schedule[1].inbounds == [schedule[3], schedule[6]] 

    # Marginal rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(x, 1.0)
    GaussianMeanPrecision(y, 0.0, 1.0)
    rf = RecognitionFactorization([x,y], ids=[:XY])
    schedule = variationalSchedule(rf.recognition_factors[:XY])
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(schedule)
    marginal_schedule = marginalSchedule(rf.recognition_factors[:XY], schedule)
    target_to_marginal_entry = ForneyLab.targetToMarginalEntry(marginal_schedule)
    assembleMarginalSchedule!(marginal_schedule, interface_to_schedule_entry, target_to_marginal_entry)
    @test marginal_schedule[3].marginal_update_rule == ForneyLab.MGaussianMeanPrecisionGGD
    @test marginal_schedule[3].marginal_id == :y_x
    @test marginal_schedule[3].inbounds == [schedule[3], schedule[1], g.nodes[:clamp_3]]
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
    schedule = sumProductSchedule(x)
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(schedule)
    assembleSchedule!(schedule, interface_to_schedule_entry, Dict())
    @test schedule[3].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test schedule[6].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
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
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    assembleSchedule!(rf.schedule, interface_to_schedule_entry, Dict())
    assembleInitialization!(rf, interface_to_schedule_entry)
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
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    assembleSchedule!(rf.schedule, interface_to_schedule_entry, Dict())
    assembleInitialization!(rf, interface_to_schedule_entry)
    @test rf.schedule[7].message_update_rule == ForneyLab.SPNonlinearIn1GG
    @test rf.schedule[3].initialize

    # Optimize
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    placeholder(x, :x)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(rf.schedule)
    assembleSchedule!(rf.schedule, interface_to_schedule_entry, Dict())
    assembleInitialization!(rf, interface_to_schedule_entry)
    @test rf.schedule[3].initialize
end

@testset "assembleAlgorithm!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    schedule = sumProductSchedule(x)
    rf.schedule = condense(flatten(schedule))
    rf.marginal_schedule = marginalSchedule(x)
    assembleAlgorithm!(rf)
    @test rf.schedule[1].schedule_index == 1
    @test rf.schedule[1].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test rf.schedule[1].inbounds == [nothing, g.nodes[:clamp_1], g.nodes[:clamp_2]]
    @test rf.schedule[2].schedule_index == 2
    @test rf.schedule[2].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test rf.schedule[2].inbounds == [nothing, g.nodes[:clamp_3], g.nodes[:clamp_4]]
    @test rf.marginal_schedule[1].marginal_id == :x
    @test rf.marginal_schedule[1].marginal_update_rule == ForneyLab.Product
    @test rf.marginal_schedule[1].inbounds == [schedule[3], schedule[6]]
end

end # module