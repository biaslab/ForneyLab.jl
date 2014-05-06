facts("GainAdditionCompositeNode") do
    context("GainAdditionCompositeNode() should initialize a GainAdditionCompositeNode with 3 interfaces") do
    end

    context("GainAdditionCompositeNode() should define an internal Equality and FixedGain node") do
    end

    context("GainAdditionCompositeNode() should point its own interfaces to the internal node interfaces") do
    end

    context("A GainAdditionCompositeNode should be able to pass a Gaussian message through its internals") do
    end

    context("A GainAdditionCompositeNode should pass a Gaussian message using custom update rules for message passing") do
        # The following tests on the update rules correspond to node 6 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        node = GainAdditionCompositeNode(2.0*eye(2)) # Defaults flag to true, telling updateNodeMessage! to use the shortcut update rules.
        context("GaussianMessage with (xi,W) parametrization") do
        end

        context("GaussianMessage with (m,W) parametrization") do
        end

        context("GaussianMessage with (m,V) parametrization") do
        end
    end
end