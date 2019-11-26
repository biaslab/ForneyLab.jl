module JuliaMessagePassingTest

using Test
using ForneyLab
import ForneyLab: assembleAlgorithm, assembleSchedule, assembleInitialization!, assembleBreaker!, assembleMarginalSchedule, assembleMessageInbound, assembleMarginalInbound

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

@testset "assembleMarginalSchedule" begin
    # Nothing rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    marginal_schedule = marginalSchedule(x)
    marginal_schedule_vect = assembleMarginalSchedule(marginal_schedule, interface_to_msg_idx)
    @test marginal_schedule_vect[1] == Dict(:marginal_id => :x,
                                            :marginal_update_rule => Nothing,
                                            :inbounds => Dict{Symbol,Any}[Dict(:schedule_index => 3)])
    # Product rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    marginal_schedule = marginalSchedule(x)
    marginal_schedule_vect = assembleMarginalSchedule(marginal_schedule, interface_to_msg_idx)
    @test marginal_schedule_vect[1] == Dict(:marginal_id => :x,
                                            :marginal_update_rule => ForneyLab.Product,
                                            :inbounds => Dict{Symbol,Any}[Dict(:schedule_index => 3), 
                                                                          Dict(:schedule_index => 6)])
    # Marginal rule
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(x, 1.0)
    GaussianMeanPrecision(y, 0.0, 1.0)
    rf = RecognitionFactorization([x,y], ids=[:XY])
    schedule = variationalSchedule(rf.recognition_factors[:XY])
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    marginal_schedule = marginalSchedule(rf.recognition_factors[:XY], schedule)
    marginal_schedule_vect = assembleMarginalSchedule(marginal_schedule, interface_to_msg_idx)
    @test marginal_schedule_vect[3] == Dict(:marginal_id => :y_x,
                                            :marginal_update_rule => ForneyLab.MGaussianMeanPrecisionGGD,
                                            :inbounds => Dict{Symbol,Any}[Dict(:schedule_index => 3), 
                                                                          Dict(:schedule_index => 1), 
                                                                          Dict(:variate_type => Univariate,
                                                                               :value => 1.0,
                                                                               :dist_or_msg => ProbabilityDistribution)])
end

@testset "assembleBreaker!" begin
    @test assembleBreaker!(Dict(), GaussianMeanPrecision, ()) == Dict(:family         => GaussianMeanPrecision,
                                                                      :initialize     => true,
                                                                      :dimensionality => ())
    @test assembleBreaker!(Dict(), Union{Gamma, Wishart}, ()) == Dict(:family         => Gamma,
                                                                      :initialize     => true,
                                                                      :dimensionality => ())
    @test assembleBreaker!(Dict(), Union{Gamma, Wishart}, (1,1)) == Dict(:family         => Wishart,
                                                                         :initialize     => true,
                                                                         :dimensionality => (1,1))
end

@testset "assembleSchedule" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    schedule_vect = assembleSchedule(schedule, interface_to_msg_idx)
    @test schedule_vect[3][:message_update_rule] == ForneyLab.SPGaussianMeanPrecisionOutNPP
    @test schedule_vect[6][:message_update_rule] == ForneyLab.SPGaussianMeanPrecisionOutNPP
end

@testset "assembleInitialization!" begin
    # Expectation propagation
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ Probit(x)
    placeholder(y, :y)
    schedule = expectationPropagationSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    rf_dict = Dict{Symbol, Any}(:id => Symbol("")) # Preset empty id
    rf_dict[:schedule] = assembleSchedule(schedule, interface_to_msg_idx)
    assembleInitialization!(rf_dict, schedule, interface_to_msg_idx)
    @test rf_dict[:schedule][5][:message_update_rule] == ForneyLab.EPProbitIn1GP
    @test rf_dict[:schedule][3][:initialize]

    # Nonlinear
    f(z) = z
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ Nonlinear(x, f)
    GaussianMeanPrecision(y, 0.0, 1.0)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    rf_dict = Dict{Symbol, Any}(:id => Symbol("")) # Preset empty id
    rf_dict[:schedule] = assembleSchedule(schedule, interface_to_msg_idx)
    assembleInitialization!(rf_dict, schedule, interface_to_msg_idx)
    @test rf_dict[:schedule][7][:message_update_rule] == ForneyLab.SPNonlinearIn1GG
    @test rf_dict[:schedule][3][:initialize]

    # Optimize
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    placeholder(x, :x)
    schedule = sumProductSchedule(x)
    interface_to_msg_idx = ForneyLab.interfaceToScheduleEntryIdx(schedule)
    rf_dict = Dict{Symbol, Any}(:id => Symbol("")) # Preset empty id
    rf_dict[:schedule] = assembleSchedule(schedule, interface_to_msg_idx)
    assembleInitialization!(rf_dict, schedule, interface_to_msg_idx)
    @test rf_dict[:schedule][3][:initialize]
end

@testset "assembleAlgorithm" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    GaussianMeanPrecision(x, 0.0, 1.0)
    schedule = sumProductSchedule(x)
    marginal_schedule = marginalSchedule(x)
    algo_dict = assembleAlgorithm(schedule, marginal_schedule)
    @test algo_dict[:schedule] == [Dict(:schedule_index => 1,
                                        :message_update_rule => ForneyLab.SPGaussianMeanPrecisionOutNPP,
                                        :inbounds => Dict{Symbol,Any}[Dict(:nothing => true), 
                                                                      Dict(:variate_type => Univariate,
                                                                           :value => 0.0,
                                                                           :dist_or_msg => Message), 
                                                                      Dict(:variate_type => Univariate,
                                                                           :value => 1.0,
                                                                           :dist_or_msg => Message)]),
                                   Dict(:schedule_index => 2,
                                        :message_update_rule => ForneyLab.SPGaussianMeanPrecisionOutNPP,
                                        :inbounds => Dict{Symbol,Any}[Dict(:nothing => true), 
                                                                      Dict(:variate_type => Univariate,
                                                                           :value => 0.0,
                                                                           :dist_or_msg => Message), 
                                                                      Dict(:variate_type => Univariate,
                                                                           :value => 1.0,
                                                                           :dist_or_msg => Message)])]
    @test algo_dict[:marginal_schedule] == [Dict(:marginal_id => :x,
                                                 :marginal_update_rule => ForneyLab.Product,
                                                 :inbounds => Dict{Symbol,Any}[Dict(:schedule_index => 1), 
                                                                               Dict(:schedule_index => 2)])]
end

end # module