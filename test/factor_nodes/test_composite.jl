module CompositeTest

using Test
using ForneyLab
using ForneyLab: @composite, outboundType, isApplicable, messagePassingSchedule
using ForneyLab: SPClamp, SPGaussianMeanVarianceOutNPP, Product, condense, flatten, step!, setTargets!


# Define new node type called StateTransition, with exposed variables called (y, x_prev, x):
@composite StateTransition (y, x_prev, x) begin
    @RV x ~ GaussianMeanVariance(x_prev, constant(1.0))
    @RV y ~ GaussianMeanVariance(x, constant(1.0))
end

@testset "@composite" begin
    @test StateTransition <: FactorNode
end

@testset "Composite node construction" begin
    g = FactorGraph()
    nd = StateTransition(Variable(), Variable(), Variable())

    # Node fields should be of correct types
    @test isa(nd.id, Symbol)
    @test isa(nd.interfaces, Vector{Interface})
    @test isa(nd.i, Dict)
    @test isa(nd.inner_graph, FactorGraph)
    @test isa(nd.terminals, Vector{Terminal})

    # Node constructor should automatically assign an id
    @test nd.id == :statetransition_1

    # Node constructor should assign interfaces to itself
    for iface in nd.interfaces
        @test ===(iface.node, nd)
    end

    # Node constructor should add node to graph
    @test g.nodes[:statetransition_1] == nd
end


#-------------
# Update rules
#-------------

# Define custom rule for sum-product message towards x
@sumProductRule(:node_type     => StateTransition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Message{Gaussian}, Nothing),
                :name          => SPStateTransitionX)

@testset "Custom SPStateTransitionX" begin
    @test SPStateTransitionX <: SumProductRule{StateTransition}
    @test outboundType(SPStateTransitionX) == Message{Gaussian}
    @test isApplicable(SPStateTransitionX, [Message{PointMass}, Message{Gaussian}, Nothing]) 
    @test !isApplicable(SPStateTransitionX, [Message{Gaussian}, Message{PointMass}, Nothing]) 
end

@testset "Composite node scheduling and algorithm compilation" begin
    fg = FactorGraph()

    x_prev = Variable(id=:x_prev)
    nd = GaussianMeanVariance(x_prev, constant(0.0), constant(1.0))
    x = Variable(id=:x)
    y = Variable(id=:y)
    cnd = StateTransition(placeholder(y, :y), x_prev, x)
    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, target_variables=Set{Variable}([x]))
    algo = InferenceAlgorithm(pfz)

    # Build SP schedule
    schedule = messagePassingSchedule(pf)
    @test length(schedule) == 5
    @test ScheduleEntry(nd.i[:m].partner, SPClamp{Univariate}) in schedule
    @test ScheduleEntry(nd.i[:v].partner, SPClamp{Univariate}) in schedule
    @test ScheduleEntry(nd.i[:out], SPGaussianMeanVarianceOutNPP) in schedule
    @test ScheduleEntry(cnd.i[:y].partner, SPClamp{Univariate}) in schedule
    @test ScheduleEntry(cnd.i[:x], SPStateTransitionX) in schedule
    pf.schedule = condense(flatten(schedule))

    # Build marginal schedule
    marginal_table = marginalTable(x)
    @test length(marginal_table) == 1
    @test marginal_table[1].target == x
    @test marginal_table[1].interfaces[1] == cnd.i[:x]
    @test marginal_table[1].marginal_update_rule == Nothing
    pf.marginal_table = marginal_table

    # Build SP algorithm for Julia execution
    ForneyLab.assembleInferenceAlgorithm!(algo)
    code = ForneyLab.algorithmSourceCode(algo)

    @test occursin("Array{Message}(undef, 2)", code)
    @test occursin("messages[1] = ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=0.0), Message(Univariate, PointMass, m=1.0))", code)
    @test occursin("messages[2] = ruleSPStateTransitionX(Message(Univariate, PointMass, m=data[:y]), messages[1], nothing)", code)
    @test occursin("marginals[:x] = messages[2].dist", code)
end

@testset "Composite node algorithm execution" begin
    # Implement custom rule for Julia execution
    ruleSPStateTransitionX(::Message{PointMass, Univariate}, ::Message{F, Univariate}, ::Nothing) where F<:Gaussian = Message(Univariate, GaussianMeanVariance, m=2.0, v=3.0) # Send some dummy message

    # Resulting algorithm ---
    function step!(marginals::Dict, data::Dict)

    messages = Array{Message}(undef, 2)

    messages[1] = ruleSPGaussianMeanVarianceOutNPP(nothing, Message(Univariate, PointMass, m=0.0), Message(Univariate, PointMass, m=1.0))
    messages[2] = ruleSPStateTransitionX(Message(Univariate, PointMass, m=data[:y]), messages[1], nothing)

    marginals[:x] = messages[2].dist

    end
    # ---

    marginals = Dict()
    data = Dict(:y => 1.0)
    step!(marginals, data)

    @test marginals[:x] == ProbabilityDistribution(Univariate, GaussianMeanVariance, m=2.0, v=3.0)
end

# Composite node definition without custom rule
@composite StateTransition2 (y, x_prev, x) begin
    @RV x ~ GaussianMeanVariance(x_prev, 1.0)
    @RV y ~ GaussianMeanVariance(x, 1.0)
end

@testset "Composite node algorithm compilation without custom rule" begin
    fg = FactorGraph()

    @RV x ~ GaussianMeanVariance(1.0, 0.0)    
    @RV y ~ GaussianMeanVariance(1.0, 0.0)
    @RV z ~ StateTransition2(x, y)
    placeholder(z, :z)

    pfz = PosteriorFactorization(fg)
    algo = messagePassingAlgorithm(x)
    code = algorithmSourceCode(algo)

    @test occursin("messages::Vector{Message}=Array{Message}(undef, 5)", code)
    @test !occursin("StateTransition", code)
end


end #module