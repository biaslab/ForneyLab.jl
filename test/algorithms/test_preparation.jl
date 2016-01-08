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

    FactorGraph()

    context("inferOutboundTypeAfterPostProcessing() should infer the correct outbound type after post-processing") do
        entry = ScheduleEntry(TerminalNode(GaussianDistribution()), 1, sumProduct!)
        entry.intermediate_outbound_type = GaussianDistribution
        entry.post_processing = sample
        @fact ForneyLab.inferOutboundTypeAfterPostProcessing(entry) --> DeltaDistribution{Float64}
        entry.post_processing = mean
        @fact ForneyLab.inferOutboundTypeAfterPostProcessing(entry) --> DeltaDistribution{Float64}
    end

    context("inferOutboundType!() should infer the correct outbound type for a terminal node") do
        entry = ScheduleEntry(TerminalNode(GaussianDistribution()), 1, sumProduct!)
        entry.inbound_types = [Void]
        ForneyLab.inferOutboundType!(entry, [sumProduct!])
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> GaussianDistribution
    end

    context("inferOutboundType!() should infer the correct outbound type after post_processing") do
        entry = ScheduleEntry(TerminalNode(GaussianDistribution()), 1, sumProduct!)
        entry.inbound_types = [Void]
        entry.post_processing = sample
        ForneyLab.inferOutboundType!(entry, [sumProduct!])
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> DeltaDistribution{Float64}
    end

    context("inferOutboundType!() should infer the correct outbound type for a node") do
        entry = ScheduleEntry(AdditionNode(), 3, sumProduct!)
        entry.inbound_types = [Message{GaussianDistribution}, Message{GaussianDistribution}, Void]
        ForneyLab.inferOutboundType!(entry, [sumProduct!])
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> GaussianDistribution
    end

    context("inferOutboundType!() should infer the correct outbound type for a Gaussian node with vmp") do
        entry = ScheduleEntry(GaussianNode(form=:precision), 2, vmp!)
        entry.inbound_types = [GaussianDistribution, Void, GaussianDistribution]
        ForneyLab.inferOutboundType!(entry, [sumProduct!, vmp!])
        @fact entry.intermediate_outbound_type --> GammaDistribution
        @fact entry.outbound_type --> GammaDistribution
    end

    context("inferOutboundType!() should infer the correct outbound type for a composite node") do
        node = initializeCompositeNode()
        entry = ScheduleEntry(node, 2, sumProduct!)
        entry.inbound_types = [GaussianDistribution, Void]
        ForneyLab.inferOutboundType!(entry, [sumProduct!, vmp!])
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> GaussianDistribution
    end

    context("inferOutboundType!() should infer the correct outbound type for a composite node with post-processing") do
        node = initializeCompositeNode()
        entry = ScheduleEntry(node, 2, sumProduct!)
        entry.inbound_types = [GaussianDistribution, Void]
        entry.post_processing = mean
        ForneyLab.inferOutboundType!(entry, [sumProduct!, vmp!])
        @fact entry.intermediate_outbound_type --> GaussianDistribution
        @fact entry.outbound_type --> DeltaDistribution{Float64}
    end

    context("buildExecute!() should pre-compile the execute field of the schedule entry") do
        # Without post-processing
        node = AdditionNode()
        node.i[:out].message = Message(GaussianDistribution())
        entry = ScheduleEntry(node, 3, sumProduct!)
        @fact isdefined(entry, :execute) --> false
        ForneyLab.buildExecute!(entry, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact typeof(entry.execute) --> Function

        # With post-processing
        node = AdditionNode()
        node.i[:out].message = Message(GaussianDistribution())
        entry = ScheduleEntry(node, 3, sumProduct!)
        entry.post_processing = mean
        entry.intermediate_outbound_type = GaussianDistribution
        entry.outbound_type = DeltaDistribution{Float64}
        @fact isdefined(entry, :execute) --> false
        ForneyLab.buildExecute!(entry, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact typeof(entry.execute) --> Function
    end
end