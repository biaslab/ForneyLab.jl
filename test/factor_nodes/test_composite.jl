module CompositeTest

using Base.Test
using ForneyLab
import ForneyLab: @composite, outboundType, isApplicable, SPConstant, SPGaussianMeanVariancePPV


# Define new node type called StateTransition, with exposed variables called (x_prev, x, y):
@composite StateTransition (x_prev, x, y) begin
    x ~ GaussianMeanVariance(x_prev, constant(1.0))
    y ~ GaussianMeanVariance(x, constant(1.0))
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
        @test is(iface.node, nd)
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
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPStateTransitionGVP)

@testset "Custom SPStateTransitionGVP" begin
    @test SPStateTransitionGVP <: SumProductRule{StateTransition}
    @test outboundType(SPStateTransitionGVP) == Message{Gaussian}
    @test isApplicable(SPStateTransitionGVP, [Message{Gaussian}, Void, Message{PointMass}]) 
    @test !isApplicable(SPStateTransitionGVP, [Message{PointMass}, Void, Message{Gaussian}]) 
end

@testset "Composite node scheduling and algorithm compilation" begin
    g = FactorGraph()

    x_prev = Variable(id=:x_prev)
    nd = GaussianMeanVariance(x_prev, constant(0.0), constant(1.0))
    x = Variable(id=:x)
    y = Variable(id=:y)
    cnd = StateTransition(x_prev, x, placeholder(y, :y))

    # Build SP schedule
    schedule = sumProductSchedule(x)
    @test length(schedule) == 5
    @test schedule[1] == ScheduleEntry(nd.i[:mean].partner, SPConstant)
    @test schedule[2] == ScheduleEntry(nd.i[:variance].partner, SPConstant)
    @test schedule[3] == ScheduleEntry(nd.i[:out], SPGaussianMeanVariancePPV)
    @test schedule[4] == ScheduleEntry(cnd.i[:y].partner, SPConstant)
    @test schedule[5] == ScheduleEntry(cnd.i[:x], SPStateTransitionGVP)

    # Build SP algorithm for Julia execution
    algo = ForneyLab.messagePassingAlgorithm(schedule, x)
    @test contains(algo, "messages = Array{Message}(2)")
    @test contains(algo, "messages[1] = ruleSPGaussianMeanVariancePPV(Message(PointMass, m=0.0), Message(PointMass, m=1.0), nothing)")
    @test contains(algo, "messages[2] = ruleSPStateTransitionGVP(messages[1], nothing, Message(PointMass, m=data[:y]))")
    @test contains(algo, "marginals[:x] = messages[2].dist")
end

@testset "Composite node algorithm execution" begin
    # Implement custom rule for Julia execution
    ruleSPStateTransitionGVP(::Message{Gaussian}, ::Void, ::Message{PointMass}) = Message(Gaussian, m=2.0, v=3.0) # Send some dummy message

    # Resulting algorithm ---
    function step!(marginals::Dict, data::Dict)

    messages = Array{Message}(2)

    messages[1] = ruleSPGaussianMeanVariancePPV(Message(PointMass, m=0.0), Message(PointMass, m=1.0), nothing)
    messages[2] = ruleSPStateTransitionGVP(messages[1], nothing, Message(PointMass, m=data[:y]))

    marginals[:x] = messages[2].dist

    end
    # ---

    marginals = Dict()
    data = Dict(:y => 1.0)
    step!(marginals, data)

    @test marginals[:x] == ProbabilityDistribution(Gaussian, m=2.0, v=3.0)
end

end #module