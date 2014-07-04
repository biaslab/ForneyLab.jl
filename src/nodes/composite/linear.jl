############################################
# LinearCompositeNode
############################################
# Description:
#   Linear variational node for parameter estimation on a linear function.
#
#         a_in  b_in   s_in
#            |  |       |
#       |-----------------|
#       |    |  |       | |
#       |    |  -->[N]<-- |
#       |    |      |     |
#       |    v      v     |
#  in1--|-->[a]--->[+]----|-->out
#       |                 |
#       |-----------------|
#
#
#
# Interface ids, (names) and supported message types:
#   1. in1:
#       GaussianMessage
#   2. a_in:
#       GaussianMessage
#   3. b_in:
#       GaussianMessage
#   4. s_in:
#       InvertedGammaMessage
#   3. out:
#       GaussianMessage
#
############################################

export LinearCompositeNode

type LinearCompositeNode <: CompositeNode
    # Basic node properties.
    # use_composite_update_rules is a flag for indicating use of shortcut update rules.
    # The linear composite node only has explicit variational update rules and thus does not need an explicit internal structure.

    use_composite_update_rules::Bool
    variational::Bool
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Helper fields filled by constructor
    in1::Interface
    a_in::Interface
    b_in::Interface
    s_in::Interface
    out::Interface

    function LinearCompositeNode(use_composite_update_rules::Bool=true, variational::Bool=true; args...)
        (name = getArgumentValue(args, :name))!=false || (name = "unnamed")

        if use_composite_update_rules == false # Check
            error("LinearCompositeNode $(name) does not support explicit internal message passing")
        end
        if variational == false # Check
            error("LinearCompositeNode $(name) only supports variational message passing")
        end

        self = new(use_composite_update_rules, variational, name, Array(Interface, 5))

        # Define the internals of the composite node
        # Initialize the composite node interfaces belonging to the composite node itself.
        for i = 1:5
            self.interfaces[i] = Interface(self)
        end
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.a_in = self.interfaces[2]
        self.b_in = self.interfaces[3]
        self.s_in = self.interfaces[4]
        self.out = self.interfaces[5]

        return self
    end
end

############################################
# Variational update functions
############################################

function updateNodeMessage!(outbound_interface_id::Int,
                            node::LinearCompositeNode,
                            inbound_messages_types::Type{Union(GaussianMessage, InverseGammaMessage)})
    # Variational update function, takes the marginals as input.
    # Derivation for the update rule can be found in the derivations notebook.

    # Ensure right parameterization
    for i = [1, 2, 3, 5]
        ensureMVParametrization!(node.interfaces[i].edge.marginal)
    end
    # Get the variables beforehand for more readable update equations
    mu_a = node.a_in.edge.marginal.m[1]
    s_a = node.a_in.edge.marginal.V[1, 1] # sigma_a^2
    mu_b = node.b_in.edge.marginal.m[1]
    s_b = node.b_in.edge.marginal.V[1, 1] # sigma_b^2
    mu_x2 = node.out.edge.marginal.m[1]
    s_x2 = node.out.edge.marginal.V[1, 1] # sigma_x2^2
    mu_x1 = node.in1.edge.marginal.m[1]
    s_x1 = node.in1.edge.marginal.V[1, 1] # sigma_x1^2
    a_s = node.s_in.edge.marginal.a
    b_s = node.s_in.edge.marginal.b

    if outbound_interface_id == 1 # in1
        nu_out = GaussianMessage( m = [(mu_a*(mu_b - mu_x2))/(s_a + mu_a^2)], V = [(a_s + 1)/(b_s*(mu_a^2 + s_a))] )
    elseif outbound_interface_id == 2 # a
        nu_out = GaussianMessage( m = [(mu_x1*(mu_b - mu_x2))/(s_x1 + mu_x1^2)], V = [(a_s + 1)/(b_s*(mu_x1^2 + s_x1))] )
    elseif outbound_interface_id == 3 # b
        nu_out = GaussianMessage( m = [mu_a*mu_x1 - mu_x2], V = [(a_s + 1)/b_s] )
    elseif outbound_interface_id == 5 # out
        nu_out = GaussianMessage( m = [mu_a*mu_x1 + mu_b], V = [(a_s + 1)/b_s] )
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = nu_out
end

function updateNodeMessage!(outbound_interface_id::Int,
                            node::LinearCompositeNode,
                            inbound_messages_types::Type{GaussianMessage})
    # Variational update function, takes the marginals as input.
    # Derivation for the update rule can be found in the derivations notebook.

    # Get the variables beforehand for more readable update equations
    mu_a = node.a_in.edge.marginal.m[1]
    s_a = node.a_in.edge.marginal.V[1, 1] # sigma_a^2
    mu_b = node.b_in.edge.marginal.m[1]
    s_b = node.b_in.edge.marginal.V[1, 1] # sigma_b^2
    mu_x2 = node.out.edge.marginal.m[1]
    s_x2 = node.out.edge.marginal.V[1, 1] # sigma_x2^2
    mu_x1 = node.in1.edge.marginal.m[1]
    s_x1 = node.in1.edge.marginal.V[1, 1] # sigma_x1^2

    if outbound_interface_id == 4 # s_N
        nu_out = InverseGammaMessage( a = -0.5, b = 0.5*((mu_x2 - mu_a*mu_x1 - mu_b)^2 - (mu_a^2)*(mu_x1^2) + (mu_a^2 + s_a)*(mu_x1^2 + s_x1) + s_b + s_x2) )
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = nu_out
end
