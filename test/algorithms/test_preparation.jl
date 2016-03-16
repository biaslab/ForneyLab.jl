#####################
# Unit tests
#####################

facts("Shared preparation methods for inference algorithms") do
    context("injectParameters!() should fill the parameters of a destination distribution with the parameters of a source distribution") do
        source = GaussianDistribution(m=4.0, V=5.0)
        destination = GaussianDistribution()
        ForneyLab.injectParameters!(destination, source)
        @fact destination.m --> 4.0
        @fact destination.V --> 5.0
        @fact is(destination, source) --> false
        @fact_throws ForneyLab.injectParameters!(destination, DeltaDistribution(4.0))
    end

    context("extractParameters() should return a dictionary of relevant type parameters") do
        call_signature_add = [AdditionNode, Type{Val{2}}, Any, Message{MvGaussianDistribution{2}}, Void, Message{MvGaussianDistribution{2}}]
        rule_signature_add = methods(sumProductRule!, call_signature_add)[1].sig.types
        @fact ForneyLab.parameters(rule_signature_add[4])[1] --> TypeVar(:dims, Union{}, Any, false)
        @fact ForneyLab.parameters(call_signature_add[4])[1] --> 2
        @fact ForneyLab.extractParameters(rule_signature_add, call_signature_add) --> Dict(TypeVar(:dims, Union{}, Any, false) => 2)

        call_signature_add_delta = [AdditionNode, Type{Val{2}}, Any, Message{MvDeltaDistribution{Float64, 2}}, Void, Message{MvDeltaDistribution{Float64, 2}}]
        rule_signature_add_delta = methods(sumProductRule!, call_signature_add_delta)[1].sig.types
        @fact ForneyLab.parameters(rule_signature_add_delta[4])[1] --> Float64
        @fact ForneyLab.parameters(call_signature_add_delta[4])[1] --> Float64
        @fact ForneyLab.parameters(rule_signature_add_delta[4])[2] --> TypeVar(:dims, Union{}, Any, false)
        @fact ForneyLab.parameters(call_signature_add_delta[4])[2] --> 2
        @fact ForneyLab.extractParameters(rule_signature_add_delta, call_signature_add_delta) --> Dict(Float64 => Float64, TypeVar(:dims, Union{}, Any, false) => 2)

        call_signature_gauss = [GaussianNode, Type{Val{2}}, Any, MvGaussianDistribution{2}, Void, MvGaussianDistribution{2}]
        rule_signature_gauss = methods(variationalRule!, call_signature_gauss)[1].sig.types
        @fact ForneyLab.parameters(rule_signature_gauss[4])[1] --> TypeVar(:dims, Union{}, Any, true)
        @fact ForneyLab.parameters(call_signature_gauss[4])[1] --> 2
        @fact ForneyLab.extractParameters(rule_signature_gauss, call_signature_gauss) --> Dict(TypeVar(:dims, Union{}, Any, true) => 2)
    end

    context("collectAllOutboundTypes() should return outbound types of applicable update rules") do
        FactorGraph()

        call_signature = [TerminalNode, Type{Val{1}}, Any, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, TerminalNode()) --> [GaussianDistribution]

        call_signature = [TerminalNode, Type{Val{1}}, Any, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, TerminalNode([1.0, 1.0])) --> [MvDeltaDistribution{Float64, 2}]

        call_signature = [EqualityNode, Type{Val{2}}, Any, Message{GaussianDistribution}, Void, Message{GaussianDistribution}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> [GaussianDistribution]

        call_signature = [EqualityNode, Type{Val{3}}, Any, Message{MvGaussianDistribution{2}}, Message{MvDeltaDistribution{Float64, 2}}, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> [MvDeltaDistribution{Float64, 2}]

        call_signature = [AdditionNode, Type{Val{2}}, Any, Message{MvGaussianDistribution{2}}, Void, Message{MvGaussianDistribution{2}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, AdditionNode()) --> [MvGaussianDistribution{2}]

        call_signature = [AdditionNode, Type{Val{2}}, Any, Message{MvDeltaDistribution{Float64, 2}}, Void, Message{MvDeltaDistribution{Float64, 2}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, AdditionNode()) --> [MvDeltaDistribution{Float64, 2}]

        call_signature = [GaussianNode{Val{:mean}, Val{:precision}}, Type{Val{2}}, Any, MvGaussianDistribution{2}, Void, MvGaussianDistribution{2}]
        @fact ForneyLab.collectAllOutboundTypes(variationalRule!, call_signature, GaussianNode(form=:precision)) --> [WishartDistribution{2}]

        call_signature = [GaussianNode{Val{:mean}, Val{:precision}}, Type{Val{3}}, Any, NormalGammaDistribution, NormalGammaDistribution, Void]
        @fact ForneyLab.collectAllOutboundTypes(variationalRule!, call_signature, GaussianNode(form=:precision)) --> [GaussianDistribution]

        call_signature = [GaussianNode{Val{:mean}, Val{:precision}}, Type{Val{1}}, Any, Void, Message{GammaDistribution}, GaussianDistribution]
        @fact ForneyLab.collectAllOutboundTypes(variationalRule!, call_signature, GaussianNode(form=:precision)) --> [StudentsTDistribution]

        call_signature = [SigmoidNode, Type{Val{1}}, Any, Message{GaussianDistribution}, Message{DeltaDistribution{Bool}}]
        @fact ForneyLab.collectAllOutboundTypes(expectationRule!, call_signature, SigmoidNode()) --> [GaussianDistribution]
    end

    context("collectAllOutboundTypes() should return outbound types for nodes with a fixed gain under inputs of different dimensions") do
        FactorGraph()

        call_signature = [GainNode, Type{Val{1}}, Any, Void, Message{MvGaussianDistribution{3}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainNode(gain=rand(3,2))) --> [MvGaussianDistribution{2}]

        call_signature = [GainNode, Type{Val{2}}, Any, Message{MvGaussianDistribution{2}}, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainNode(gain=rand(3,2))) --> [MvGaussianDistribution{3}]

        call_signature = [GainAdditionNode, Type{Val{1}}, Any, Void, Message{MvGaussianDistribution{3}}, Message{MvGaussianDistribution{3}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainAdditionNode(rand(3,2))) --> [MvGaussianDistribution{2}]

        call_signature = [GainEqualityNode, Type{Val{3}}, Any, Message{MvGaussianDistribution{2}}, Message{MvGaussianDistribution{2}}, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainEqualityNode(rand(3,2))) --> [MvGaussianDistribution{3}]
    end

    FactorGraph()

    context("inferOutboundTypeAfterPostProcessing() should infer the correct outbound type after post-processing") do
        entry = ScheduleEntry(TerminalNode(GaussianDistribution()), 1, sumProductRule!)
        entry.intermediate_outbound_type = GaussianDistribution
        entry.post_processing = sample
        @fact ForneyLab.inferOutboundTypeAfterPostProcessing(entry) --> DeltaDistribution{Float64}
        entry.post_processing = mean
        @fact ForneyLab.inferOutboundTypeAfterPostProcessing(entry) --> DeltaDistribution{Float64}
    end

    context("inferOutboundType!() should infer the correct outbound type for a terminal node") do
        entry = ScheduleEntry(TerminalNode(GaussianDistribution()), 1, sumProductRule!)
        entry.inbound_types = [Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> GaussianDistribution
    end

    context("inferOutboundType!() should infer the correct outbound type after post_processing") do
        entry = ScheduleEntry(TerminalNode(GaussianDistribution()), 1, sumProductRule!)
        entry.inbound_types = [Void]
        entry.post_processing = sample
        ForneyLab.inferOutboundType!(entry)
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> DeltaDistribution{Float64}
    end

    context("inferOutboundType!() should infer the correct outbound type for a node") do
        entry = ScheduleEntry(AdditionNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{GaussianDistribution}, Message{GaussianDistribution}, Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> GaussianDistribution
    end

    context("inferOutboundType!() should infer the correct outbound type for a Gaussian node with vmp") do
        entry = ScheduleEntry(GaussianNode(form=:precision), 2, variationalRule!)
        entry.inbound_types = [GaussianDistribution, Void, GaussianDistribution]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.intermediate_outbound_type --> GammaDistribution
        @fact entry.outbound_type --> GammaDistribution
    end

    context("buildExecute!() should pre-compile the execute field of the schedule entry") do
        # Without post-processing
        node = AdditionNode()
        node.i[:out].message = Message(GaussianDistribution())
        entry = ScheduleEntry(node, 3, sumProductRule!)
        @fact isdefined(entry, :execute) --> false
        ForneyLab.buildExecute!(entry, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact typeof(entry.execute) --> Function

        # With post-processing
        node = AdditionNode()
        node.i[:out].message = Message(GaussianDistribution())
        entry = ScheduleEntry(node, 3, sumProductRule!)
        entry.post_processing = mean
        entry.intermediate_outbound_type = GaussianDistribution
        entry.outbound_type = DeltaDistribution{Float64}
        @fact isdefined(entry, :execute) --> false
        ForneyLab.buildExecute!(entry, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact typeof(entry.execute) --> Function
    end
end
