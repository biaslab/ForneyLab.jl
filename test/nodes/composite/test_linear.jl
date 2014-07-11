#####################
# Unit tests
#####################

facts("LinearCompositeNode unit tests") do
    context("LinearCompositeNode() should initialize a LinearCompositeNode with 5 interfaces") do
        node = LinearCompositeNode()
        @fact typeof(node) => LinearCompositeNode
        @fact length(node.interfaces) => 5
        @fact node.in1 => node.interfaces[1]
        @fact node.a_in => node.interfaces[2]
        @fact node.b_in => node.interfaces[3]
        @fact node.noise_in => node.interfaces[4]
        @fact node.out => node.interfaces[5]
        @fact node.variational => true # default variational to true
        @fact node.use_composite_update_rules => true # default use_composite_update_rules to true
    end
end

#####################
# Integration tests
#####################

facts("LinearCompositeNode integration tests") do
    context("LinearCompositeNode should propagate a backward variational message to in1") do
        lin_node = initializeLinearCompositeNode([uninformative(GaussianDistribution), GaussianDistribution(m=[2.0], V=[0.0]), GaussianDistribution(m=[0.5], V=[0.0]), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=[2.5], V=[0.0])])
        msg = ForneyLab.updateNodeMessage!(1, lin_node, Union(GaussianDistribution, InverseGammaDistribution))
        @fact msg.value.m => [1.0]
    end

    context("LinearCompositeNode should propagate a backward variational message to a_in") do
        lin_node = initializeLinearCompositeNode([GaussianDistribution(m=[1.0], V=[0.0]), uninformative(GaussianDistribution), GaussianDistribution(m=[0.5], V=[0.0]), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=[2.5], V=[0.0])])
        msg = ForneyLab.updateNodeMessage!(2, lin_node, Union(GaussianDistribution, InverseGammaDistribution))
        @fact msg.value.m => [2.0]
    end

    context("LinearCompositeNode should propagate a backward variational message to b_in") do
        lin_node = initializeLinearCompositeNode([GaussianDistribution(m=[1.0], V=[0.0]), GaussianDistribution(m=[2.0], V=[0.0]), uninformative(GaussianDistribution), InverseGammaDistribution(a=10000.0, b=19998.0), GaussianDistribution(m=[2.5], V=[0.0])])
        msg = ForneyLab.updateNodeMessage!(3, lin_node, Union(GaussianDistribution, InverseGammaDistribution))
        @fact msg.value.m => [0.5]
    end

    context("LinearCompositeNode should propagate a backward variational message to noise_in") do
        lin_node = initializeLinearCompositeNode([GaussianDistribution(m=[1.0], V=[0.0]), GaussianDistribution(m=[2.0], V=[0.0]), GaussianDistribution(m=[0.5], V=[0.0]), uninformative(InverseGammaDistribution), GaussianDistribution(m=[2.5], V=[0.0])])
        msg = ForneyLab.updateNodeMessage!(4, lin_node, Union(GaussianDistribution))
        @fact msg.value.a => -0.5
        @fact msg.value.b => 0.0
    end

    context("LinearCompositeNode should propagate a forward variational message to out") do
        lin_node = initializeLinearCompositeNode([GaussianDistribution(m=[1.0], V=[0.0]), GaussianDistribution(m=[2.0], V=[0.0]), GaussianDistribution(m=[0.5], V=[0.0]), InverseGammaDistribution(a=10000.0, b=19998.0), uninformative(GaussianDistribution)])
        msg = ForneyLab.updateNodeMessage!(5, lin_node, Union(GaussianDistribution, InverseGammaDistribution))
        @fact msg.value.m => [2.5]
    end
end