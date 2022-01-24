module MessagePassingAssemblersTest

using Test
using ForneyLab
using ForneyLab: assembleBreaker!, assembleClamp!, assembleInferenceAlgorithm!, assemblePosteriorFactor!, assembleSchedule!, assembleInitialization!, assembleMarginalTable!, condense, flatten, setTargets!, messagePassingSchedule

@testset "assembleClamp!" begin
    g = FactorGraph()
    nd = Clamp(Variable(), 1.0)
    assembleClamp!(nd, Distribution)    
    @test nd.dist_or_msg == Distribution
end

@testset "assembleBreaker!" begin
    breaker_entry = ScheduleEntry()
    assembleBreaker!(breaker_entry, Gaussian{Precision}, ())
    @test breaker_entry.family == Gaussian{Precision}
    @test breaker_entry.initialize == true
    @test breaker_entry.dimensionality == ()

    breaker_entry = ScheduleEntry()
    assembleBreaker!(breaker_entry, Union{Gamma, Wishart}, (1,1))
    @test breaker_entry.family == Union{Gamma, Wishart}
    @test breaker_entry.initialize == true
    @test breaker_entry.dimensionality == (1,1)
end

@testset "assembleSchedule!" begin
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    Gaussian{Precision}(x, 0.0, 1.0)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    pf.schedule = messagePassingSchedule(pf)
    algo = InferenceAlgorithm(pfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(pf)
    @test pf.schedule[3].message_update_rule == ForneyLab.SPGaussianPrecisionOutNPP
    @test pf.schedule[6].message_update_rule == ForneyLab.SPGaussianPrecisionOutNPP
end

@testset "assembleInitialization!" begin
    # Expectation propagation
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    @RV y ~ Probit(x)
    placeholder(y, :y)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    pf.ep_sites = Set{Interface}([fg.nodes[:probit_1].i[:in1]])
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    pf.schedule = messagePassingSchedule(pf)
    algo = InferenceAlgorithm(pfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(pf)
    assembleInitialization!(pf)
    @test pf.schedule[5].message_update_rule == ForneyLab.EPProbitIn1PG
    @test pf.schedule[3].initialize

    # Delta
    f(z) = z
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    @RV y ~ Delta{Unscented}(x, g=f)
    Gaussian{Precision}(y, 0.0, 1.0)
    
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    pf.schedule = messagePassingSchedule(pf)
    algo = InferenceAlgorithm(pfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(pf)
    assembleInitialization!(pf)
    @test pf.schedule[7].message_update_rule == ForneyLab.SPDeltaUTIn1GG
    @test pf.schedule[3].initialize

    # Optimize
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    placeholder(x, :x)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    pf.schedule = messagePassingSchedule(pf)
    algo = InferenceAlgorithm(pfz)
    algo.target_to_marginal_entry = Dict()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    assembleSchedule!(pf)
    assembleInitialization!(pf)
    @test pf.schedule[3].initialize
end

@testset "assembleMarginalTable!" begin
    # Nothing rule
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    pf.schedule = messagePassingSchedule(pf)
    pf.marginal_table = marginalTable(x)
    algo = InferenceAlgorithm(pfz)
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleMarginalTable!(pf)
    @test pf.marginal_table[1].marginal_update_rule == Nothing
    @test pf.marginal_table[1].marginal_id == :x
    @test pf.marginal_table[1].inbounds == [pf.schedule[3]]

    # Product rule
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    Gaussian{Precision}(x, 0.0, 1.0)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    pf.schedule = messagePassingSchedule(pf)
    pf.marginal_table = marginalTable(x)
    algo = InferenceAlgorithm(pfz)
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleMarginalTable!(pf)
    @test pf.marginal_table[1].marginal_update_rule == ForneyLab.Product
    @test pf.marginal_table[1].marginal_id == :x
    @test pf.marginal_table[1].inbounds == [pf.schedule[3], pf.schedule[6]] 

    # Marginal rule
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    @RV y ~ Gaussian{Precision}(x, 1.0)
    Gaussian{Precision}(y, 0.0, 1.0)
    pfz = PosteriorFactorization([x,y], ids=[:XY])
    pf = pfz.posterior_factors[:XY]
    setTargets!(pf, pfz, external_targets=true)
    pf.schedule = messagePassingSchedule(pf)
    pf.marginal_table = marginalTable(pf)
    algo = InferenceAlgorithm(pfz)
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assembleMarginalTable!(pf)
    @test pf.marginal_table[3].marginal_update_rule == ForneyLab.MGaussian{Precision}GGD
    @test pf.marginal_table[3].marginal_id == :y_x
    @test length(pf.marginal_table[3].inbounds) == 3
end

@testset "assemblePosteriorFactor!" begin
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    Gaussian{Precision}(x, 0.0, 1.0)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    schedule = messagePassingSchedule(pf)
    pf.schedule = condense(flatten(schedule))
    pf.marginal_table = marginalTable(x)
    algo = InferenceAlgorithm()
    algo.interface_to_schedule_entry = ForneyLab.interfaceToScheduleEntry(algo)
    algo.target_to_marginal_entry = ForneyLab.targetToMarginalEntry(algo)
    assemblePosteriorFactor!(pf)
    @test pf.schedule[1].schedule_index == 1
    @test pf.schedule[1].message_update_rule == ForneyLab.SPGaussianPrecisionOutNPP
    @test pf.schedule[1].inbounds == [nothing, fg.nodes[:clamp_1], fg.nodes[:clamp_2]]
    @test pf.schedule[2].schedule_index == 2
    @test pf.schedule[2].message_update_rule == ForneyLab.SPGaussianPrecisionOutNPP
    @test pf.schedule[2].inbounds == [nothing, fg.nodes[:clamp_3], fg.nodes[:clamp_4]]
    @test pf.marginal_table[1].marginal_id == :x
    @test pf.marginal_table[1].marginal_update_rule == ForneyLab.Product
    @test pf.marginal_table[1].inbounds == [schedule[3], schedule[6]]
end

@testset "assembleInferenceAlgorithm!" begin
    fg = FactorGraph()
    @RV x ~ Gaussian{Precision}(0.0, 1.0)
    Gaussian{Precision}(x, 0.0, 1.0)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    schedule = messagePassingSchedule(pf)
    pf.schedule = condense(flatten(schedule))
    pf.marginal_table = marginalTable(x)
    algo = InferenceAlgorithm(pfz)
    assembleInferenceAlgorithm!(algo)

    @test length(algo.interface_to_schedule_entry) == 2
    @test length(algo.target_to_marginal_entry) == 1
end

end # module