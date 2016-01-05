#####################
# Unit tests
#####################

facts("SigmoidNode unit tests") do
    context("Construction") do
        FactorGraph()
        sig_node = SigmoidNode(id=:node)
        @fact length(sig_node.interfaces) --> 2
        @fact sig_node.i[:real] --> sig_node.interfaces[1]
        @fact sig_node.i[:bin] --> sig_node.interfaces[2]
    end

    context("sumProduct! rules") do
        # Forward message, DeltaDistribution on i[:real]
        validateOutboundMessage(SigmoidNode(),
                                2,
                                [Message(DeltaDistribution(0.0)), nothing],
                                BernoulliDistribution(0.5))
        validateOutboundMessage(SigmoidNode(),
                                2,
                                [Message(DeltaDistribution(1.0)), nothing],
                                BernoulliDistribution(ForneyLab.Φ(1.0)))
        # Forward message, GaussianDistribution on i[:real]
        validateOutboundMessage(SigmoidNode(),
                                2,
                                [Message(GaussianDistribution(m=1.0, V=0.5)), nothing],
                                BernoulliDistribution(ForneyLab.Φ(1/sqrt(1+0.5))))
    end

    context("ep! rules") do
        # TODO: update
        @fact true --> false

        # FactorGraph()
        # sig_node = SigmoidNode(id=:sigmoid)
        # cavity = GaussianDistribution(m=0.0, V=1.0)
        # Edge(TerminalNode(cavity, id=:t_real), sig_node.i[:real])
        # Edge(sig_node.i[:bin], TerminalNode(DeltaDistribution(true), id=:t_bin))

        # # Backward message, DeltaDistribution on i[:bin]
        # n(:t_bin).value = DeltaDistribution(true)
        # @fact ForneyLab.ep!(sig_node, 1, Message(n(:t_real).value), Message(n(:t_bin).value)) --> (:sigmoid_backward_gaussian_expectation, sig_node.i[:real].message)
        # ensureParameters!(sig_node.i[:real].message.payload, (:m, :V))
        # @fact sig_node.i[:real].message.payload.m --> roughly(1.7724, atol=1e-4)
        # @fact sig_node.i[:real].message.payload.V --> roughly(2.1415, atol=1e-4)

        # n(:t_bin).value = DeltaDistribution(false)
        # ForneyLab.ep!(sig_node, 1, Message(n(:t_real).value), Message(n(:t_bin).value))
        # ensureParameters!(sig_node.i[:real].message.payload, (:m, :V))
        # @fact sig_node.i[:real].message.payload.m --> roughly(-1.7724, atol=1e-4)
        # @fact sig_node.i[:real].message.payload.V --> roughly(2.1415, atol=1e-4)
        # ensureParameters!(sig_node.i[:real].edge.marginal, (:m, :V))
        # mean_marg_hard = sig_node.i[:real].edge.marginal.m

        # # Backward message, BernoulliDistribution on i[:bin]
        # # Uninformative data should result in marginal (almost) equal to cavity distribution and vague backward message
        # n(:t_bin).value = BernoulliDistribution(0.5) # uninformative data
        # ForneyLab.ep!(sig_node, 1, Message(n(:t_real).value), Message(n(:t_bin).value))
        # ensureParameters!(sig_node.i[:real].message.payload, (:m, :V))
        # ensureParameters!(sig_node.i[:real].edge.marginal, (:m, :V))
        # @fact sig_node.i[:real].message.payload.m --> roughly(0.0, atol=1e-4)
        # @fact sig_node.i[:real].message.payload.V --> greater_than(huge/2)
        # @fact sig_node.i[:real].message.payload.V --> less_than(Inf)
        # @fact sig_node.i[:real].edge.marginal.m --> roughly(n(:t_real).value.m)
        # @fact sig_node.i[:real].edge.marginal.V --> roughly(n(:t_real).value.V)

        # # Backward message, BernoulliDistribution on i[:bin]
        # n(:t_bin).value = BernoulliDistribution(0.2) # softbit
        # ForneyLab.ep!(sig_node, 1, Message(n(:t_real).value), Message(n(:t_bin).value))
        # ensureParameters!(sig_node.i[:real].message.payload, (:m, :V))
        # ensureParameters!(sig_node.i[:real].edge.marginal, (:m, :V))
        # @fact sig_node.i[:real].message.payload.m --> roughly(-2.9540, atol=1e-4)
        # @fact sig_node.i[:real].message.payload.V --> roughly(7.7266, atol=1e-4)
        # @fact sig_node.i[:real].edge.marginal.m --> greater_than(mean_marg_hard) # softbit should result in posterior closer to prior compared to 'hard' bit
    end
end
