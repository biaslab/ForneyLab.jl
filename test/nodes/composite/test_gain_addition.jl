facts("GainAdditionCompositeNode") do
    context("GainAdditionCompositeNode() should initialize a GainAdditionCompositeNode with 3 interfaces") do
        node = GainAdditionCompositeNode()
        @fact typeof(node) => GainAdditionCompositeNode
        @fact length(node.interfaces) => 3
        @fact node.in1 => node.interfaces[1]
        @fact node.in2 => node.interfaces[2]
        @fact node.out => node.interfaces[3]
        @fact typeof(node.A) => Array{Float64, 2}
    end

    context("GainAdditionCompositeNode() should define an internal AdditionNode and FixedGainNode") do
        node = GainAdditionCompositeNode([5.0], false)
        @fact typeof(node.addition_node) => AdditionNode
        @fact typeof(node.fixed_gain_node) => FixedGainNode
        @fact node.fixed_gain_node.A => reshape([5.0], 1, 1)
    end

    context("GainAdditionCompositeNode() should point its own interfaces to the internal node interfaces") do
        node = GainAdditionCompositeNode([1.0], false)
        @fact node.in1 => node.fixed_gain_node.interfaces[1]
        @fact node.in2 => node.addition_node.interfaces[2]
        @fact node.out => node.addition_node.out
    end

    context("GainAdditionCompositeNode should be able to pass GaussianMessages: using shortcut rules or internal graph should yield same result (m,V) parametrization") do
        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->
        #        |_______|
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)
        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], V=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(m=[1.0, 2.0], V=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.in2)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.out)
        @fact is(msg_shortcut, node_gain_addition_composite.out.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.in2)
        msg_internal = calculateMessage!(node_gain_addition_internal.out)
        @fact is(msg_internal, node_gain_addition_internal.out.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], V=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(m=[1.0, 2.0], V=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in2)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in1)
        @fact is(msg_shortcut, node_gain_addition_composite.in1.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in2)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in1)
        @fact is(msg_internal, node_gain_addition_internal.in1.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #   -----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], V=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(m=[1.0, 2.0], V=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in2)
        @fact is(msg_shortcut, node_gain_addition_composite.in2.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in2)
        @fact is(msg_internal, node_gain_addition_internal.in2.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true
    end

    context("GainAdditionCompositeNode should be able to pass GaussianMessages: using shortcut rules or internal graph should yield same result (m,W) parametrization") do
        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->
        #        |_______|
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)
        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], W=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(m=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.in2)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.out)
        @fact is(msg_shortcut, node_gain_addition_composite.out.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.in2)
        msg_internal = calculateMessage!(node_gain_addition_internal.out)
        @fact is(msg_internal, node_gain_addition_internal.out.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], W=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(m=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in2)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in1)
        @fact is(msg_shortcut, node_gain_addition_composite.in1.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in2)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in1)
        @fact is(msg_internal, node_gain_addition_internal.in1.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #   -----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], W=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(m=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in2)
        @fact is(msg_shortcut, node_gain_addition_composite.in2.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in2)
        @fact is(msg_internal, node_gain_addition_internal.in2.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true
    end

    context("GainAdditionCompositeNode should be able to pass GaussianMessages: using shortcut rules or internal graph should yield same result (xi,W) parametrization") do
        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->
        #        |_______|
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)
        node_c1 = ConstantNode(GaussianMessage(xi=[0.0, 0.0], W=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(xi=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.in2)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.out)
        @fact is(msg_shortcut, node_gain_addition_composite.out.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.in2)
        msg_internal = calculateMessage!(node_gain_addition_internal.out)
        @fact is(msg_internal, node_gain_addition_internal.out.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(xi=[0.0, 0.0], W=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(xi=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in2)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in1)
        @fact is(msg_shortcut, node_gain_addition_composite.in1.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in2)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in1)
        @fact is(msg_internal, node_gain_addition_internal.in1.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #   -----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(xi=[0.0, 0.0], W=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(xi=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in2)
        @fact is(msg_shortcut, node_gain_addition_composite.in2.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in2)
        @fact is(msg_internal, node_gain_addition_internal.in2.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true
    end
    context("GainAdditionCompositeNode should be able to pass GaussianMessages: using shortcut rules or internal graph should yield same result (different parametrizations)") do
        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->
        #        |_______|
        A = reshape([2.0, 3.0, 3.0, 2.0], 2, 2)
        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], V=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(xi=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.in2)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.out)
        @fact is(msg_shortcut, node_gain_addition_composite.out.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.in2)
        msg_internal = calculateMessage!(node_gain_addition_internal.out)
        @fact is(msg_internal, node_gain_addition_internal.out.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #[N]-----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], V=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(xi=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in2)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in1)
        @fact is(msg_shortcut, node_gain_addition_composite.in1.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in2)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in1)
        @fact is(msg_internal, node_gain_addition_internal.in1.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true

        #           [N]
        #            | in1
        #            |
        #        ____|____
        #        |   v   |
        #        |  [A]  |
        #        |   |   |
        #    in2 |   v   | out
        #   -----|->[+]--|---->[N]
        #        |_______|

        node_c1 = ConstantNode(GaussianMessage(m=[0.0, 0.0], V=eye(2,2)))
        node_c2 = ConstantNode(GaussianMessage(xi=[1.0, 2.0], W=2.0*eye(2,2)))

        #Shortcut rules
        node_gain_addition_composite = GainAdditionCompositeNode(A, true)
        Edge(node_c1.out, node_gain_addition_composite.in1)
        Edge(node_c2.out, node_gain_addition_composite.out)
        msg_shortcut = calculateMessage!(node_gain_addition_composite.in2)
        @fact is(msg_shortcut, node_gain_addition_composite.in2.message) => true

        #Internal graph
        node_gain_addition_internal = GainAdditionCompositeNode(A, false)
        Edge(node_c1.out, node_gain_addition_internal.in1)
        Edge(node_c2.out, node_gain_addition_internal.out)
        msg_internal = calculateMessage!(node_gain_addition_internal.in2)
        @fact is(msg_internal, node_gain_addition_internal.in2.message) => true

        ensureMVParametrization!(msg_shortcut)
        ensureMVParametrization!(msg_internal)
        @fact isApproxEqual(msg_internal.m, msg_shortcut.m) => true
        @fact isApproxEqual(msg_internal.V, msg_shortcut.V) => true
    end
end