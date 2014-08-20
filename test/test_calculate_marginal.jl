facts("Marginal calculation unit tests") do
    context("Marginal calculation for the GaussianNode") do
        node = GaussianNode(form="precision")
        m_edge = Edge(MockNode(Message(GaussianDistribution())).out, node.mean)
        gam_edge = Edge(MockNode(Message(GammaDistribution())).out, node.precision)
        calculateMarginal!(node)
        @fact node.marginal => NormalGammaDistribution(m=0.0, W=1.0, a=1.0, b=1.0) 
    end

    context("Marginal calculation for the combination of a Gaussian and student's t-distribution") do
        edge = Edge(MockNode(Message(GaussianDistribution())).out, MockNode(Message(StudentsTDistribution())).out)
        calculateMarginal!(edge)
        @fact edge.marginal => GaussianDistribution(m=0.0, W=3.0) 
    end
end