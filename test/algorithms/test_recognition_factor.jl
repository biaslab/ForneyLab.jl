module RecognitionFactorTest

using Test
using ForneyLab

import ForneyLab: nodesConnectedToExternalEdges, Cluster, hasCollider, assembleAlgorithm!, assembleSchedule!, assembleInitialization!, assembleMarginalTable!, condense, flatten

@testset "RecognitionFactor" begin
    g = FactorGraph()
    @RV m ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV w ~ Gamma(constant(1.0), constant(1.0))
    y = Variable[]
    for i = 1:3
        @RV y_i ~ GaussianMeanPrecision(m, w)
        placeholder(y_i, :y, index=i)
        push!(y, y_i)
    end

    algo = Algorithm()
    q_m = RecognitionFactor(m)
    @test q_m.id == :recognitionfactor_1
    @test q_m.variables == Set([m])
    @test q_m.clusters == Set{Cluster}()
    @test q_m.internal_edges == edges(m)
    @test algo.recognition_factors[:recognitionfactor_1] === q_m

    q_w = RecognitionFactor(w)
    @test q_w.id == :recognitionfactor_2
    @test q_w.variables == Set([w])
    @test q_w.clusters == Set{Cluster}()
    @test q_w.internal_edges == edges(w)
    @test algo.recognition_factors[:recognitionfactor_2] === q_w

    # Joint factorizations
    q_m_w = RecognitionFactor([m, w])
    @test q_m_w.id == :recognitionfactor_3
    @test q_m_w.variables == Set([m, w])
    @test length(q_m_w.clusters) == 3 
    @test q_m_w.internal_edges == edges(Set([m, w]))
    @test algo.recognition_factors[:recognitionfactor_3] === q_m_w

    q_y = RecognitionFactor(y)
    @test q_y.id == :recognitionfactor_4
    @test q_y.variables == Set(y)
    @test q_y.clusters == Set{Cluster}()
    @test q_y.internal_edges == edges(Set(y))
    @test algo.recognition_factors[:recognitionfactor_4] === q_y
end

@testset "hasCollider()" begin

    # [N]--->[+]<---[N]
    #     a   |   b
    #         v c
    #        [N]
    #         |
    #         v d
    #         ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(0.0, 1.0)
    @RV c = a + b
    @RV d ~ GaussianMeanVariance(c, 1.0)
    placeholder(d, :d)
    algo = Algorithm()
    q = RecognitionFactor([a,b,c])
    @test hasCollider(q) == true

    # [N]--->[=]<---[N]
    #     a   |   b
    #         v c
    #        [N]
    #         |
    #         v d
    #         ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(0.0, 1.0)
    @RV c = equal(a, b)
    @RV d ~ GaussianMeanVariance(c, 1.0)
    placeholder(d, :d)
    algo = Algorithm()
    q = RecognitionFactor([a,b,c])
    @test hasCollider(q) == true

    # [N]--->[=]--->[N]--->■
    #         |  a      b
    #         v 
    #        [N]
    #         |
    #         v c
    #         ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(a, 1.0)
    @RV c ~ GaussianMeanVariance(a, 1.0)
    placeholder(b, :b)
    placeholder(c, :c)
    algo = Algorithm()
    q = RecognitionFactor(a)
    @test hasCollider(q) == false

    # [N]--->[N]--->[N]--->■
    #     a      b      c

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(a, 1.0)
    @RV c ~ GaussianMeanVariance(b, 1.0)
    placeholder(c, :c)
    algo = Algorithm()
    q = RecognitionFactor([a,b])
    @test hasCollider(q) == false

    # [N]--->[N]--->■
    #     a      b 
    #
    # [N]--->[N]--->■
    #     c      d    

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(a, 1.0)
    @RV c ~ GaussianMeanVariance(0.0, 1.0)
    @RV d ~ GaussianMeanVariance(c, 1.0)
    placeholder(b, :b)
    placeholder(d, :d)
    algo = Algorithm()
    q = RecognitionFactor([a,c])
    @test hasCollider(q) == false

    # [N]--->[+]<---[N]
    #     a   |   b
    #         v c
    # [N]<---[=]--->[N]
    #  |             |
    #  v d           v e 
    #  ■             ■

    g = FactorGraph()
    @RV a ~ GaussianMeanVariance(0.0, 1.0)
    @RV b ~ GaussianMeanVariance(0.0, 1.0)
    @RV c = a + b
    @RV d ~ GaussianMeanVariance(c, 1.0)
    @RV e ~ GaussianMeanVariance(c, 1.0)
    placeholder(d, :d)
    placeholder(e, :e)
    algo = Algorithm()
    q = RecognitionFactor([a,b,c])
    @test hasCollider(q) == true
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