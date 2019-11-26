module JuliaMessagePassingTest

using Test
using ForneyLab
import ForneyLab: assembleAlgorithm!, assembleSchedule!, assembleInitialization!, assembleBreaker!, assembleMarginalSchedule!, assembleMessageInbound, assembleMarginalInbound, condense, flatten

@testset "assembleMarginalInbound" begin
    g = FactorGraph()
    nd = Clamp(Variable(), 1.0)
    inbound_dict = assembleMarginalInbound(nd)    

    @test inbound_dict == Dict(:variate_type => Univariate,
                               :value        => 1.0,
                               :dist_or_msg  => ProbabilityDistribution)
end

@testset "assembleMessageInbound" begin
    g = FactorGraph()
    nd = Clamp(Variable(), 1.0)
    inbound_dict = assembleMessageInbound(nd)    

    @test inbound_dict == Dict(:variate_type => Univariate,
                               :value        => 1.0,
                               :dist_or_msg  => Message)
end

@testset "assembleMarginalSchedule!" begin
    # Nothing rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    marginal_schedule = marginalSchedule(x)
    assembleMarginalSchedule!(marginal_schedule, interface_to_msg_idx)
    @test marginal_schedule[1].marginal_update_rule == Nothing
    @test marginal_schedule[1].marginal_id == :x
    @test marginal_schedule[1].inbounds == Dict{Symbol,Any}[Dict(:schedule_index => 3)]

    # Product rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    marginal_schedule = marginalSchedule(x)
    assembleMarginalSchedule!(marginal_schedule, interface_to_msg_idx)
    @test marginal_schedule[1].marginal_update_rule == ForneyLab.Product
    @test marginal_schedule[1].marginal_id == :x
    @test marginal_schedule[1].inbounds == Dict{Symbol,Any}[Dict(:schedule_index => 3), 
                                                            Dict(:schedule_index => 6)]
    # Marginal rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(x, 1.0)
    GaussianMeanPrecision(y, 0.0, 1.0)
    rf = RecognitionFactorization([x,y], ids=[:XY])
    schedule = variationalSchedule(rf.recognition_factors[:XY])
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    marginal_schedule = marginalSchedule(rf.recognition_factors[:XY], schedule)
    assembleMarginalSchedule!(marginal_schedule, interface_to_msg_idx)
    @test marginal_schedule[3].marginal_update_rule == ForneyLab.MGaussianMeanPrecisionGGD
    @test marginal_schedule[3].marginal_id == :y_x
    @test marginal_schedule[3].inbounds == Dict{Symbol,Any}[Dict(:schedule_index => 3), 
                                                            Dict(:schedule_index => 1), 
                                                            Dict(:variate_type => Univariate,
                                                                 :value => 1.0,
                                                                 :dist_or_msg => ProbabilityDistribution)]
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
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    assembleSchedule!(schedule, interface_to_msg_idx)
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
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(rf.schedule)
    assembleSchedule!(rf.schedule, interface_to_msg_idx)
    assembleInitialization!(rf, interface_to_msg_idx)
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
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(rf.schedule)
    assembleSchedule!(rf.schedule, interface_to_msg_idx)
    assembleInitialization!(rf, interface_to_msg_idx)
    @test rf.schedule[7].message_update_rule == ForneyLab.SPNonlinearIn1GG
    @test rf.schedule[3].initialize

    # Optimize
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    placeholder(x, :x)
    rfz = RecognitionFactorization()
    rf = RecognitionFactor(rfz)
    rf.schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(rf.schedule)
    assembleSchedule!(rf.schedule, interface_to_msg_idx)
    assembleInitialization!(rf, interface_to_msg_idx)
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
    @test rf.schedule[1].inbounds == Dict{Symbol,Any}[Dict(:nothing => true), 
                                                      Dict(:variate_type => Univariate,
                                                           :value => 0.0,
                                                           :dist_or_msg => Message), 
                                                      Dict(:variate_type => Univariate,
                                                           :value => 1.0,
                                                           :dist_or_msg => Message)]
    @test rf.schedule[2].schedule_index == 2
    @test rf.schedule[2].message_update_rule == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test rf.schedule[2].inbounds == Dict{Symbol,Any}[Dict(:nothing => true), 
                                                      Dict(:variate_type => Univariate,
                                                           :value => 0.0,
                                                           :dist_or_msg => Message), 
                                                      Dict(:variate_type => Univariate,
                                                           :value => 1.0,
                                                           :dist_or_msg => Message)]
    @test rf.marginal_schedule[1].marginal_id == :x
    @test rf.marginal_schedule[1].marginal_update_rule == ForneyLab.Product
    @test rf.marginal_schedule[1].inbounds == Dict{Symbol,Any}[Dict(:schedule_index => 1), 
                                                               Dict(:schedule_index => 2)]
end

end # module