############################################
# LinearCompositeNode
############################################
# Description:
#   Linear variational node for parameter estimation on a linear function.
#
#         slope offset  noise
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
#       GaussianDistribution (marginal)
#   2. slope:
#       GaussianDistribution (marginal)
#   3. offset:
#       GaussianDistribution (marginal)
#   4. noise:
#       InvertedGammaDistribution (marginal)
#       GammaDistribution (marginal)
#   3. out:
#       GaussianDistribution (marginal)
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
    slope::Interface
    offset::Interface
    noise::Interface
    out::Interface

    function LinearCompositeNode(use_composite_update_rules::Bool=true, variational::Bool=true; name = "unnamed", form::ASCIIString="moment", args...)
        if use_composite_update_rules == false # Check
            error("LinearCompositeNode $(name) does not support explicit internal message passing")
        end
        if variational == false # Check
            error("LinearCompositeNode $(name) only supports variational message passing")
        end

        self = new(use_composite_update_rules, variational, name, Array(Interface, 5))

        args = Dict(zip(args...)...) # Cast args to dictionary
        # Set up the interfaces
        param_list = [:in1, :slope, :offset, :noise, :out]
        for i = 1:length(param_list)
            self.interfaces[i] = Interface(self) # Construct interface
            setfield(self, param_list[i], self.interfaces[i]) # Set named interfaces

            # Clamp parameter values when given as argument
            if haskey(args, param_list[i])
                Edge(ForneyLab.ClampNode(Message(args[param_list[i]])).out, getfield(self, param_list[i]), typeof(args[param_list[i]])) # Connect clamp node
            end
        end

        return self
    end
end


############################################
# Variational update functions
############################################

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GaussianDistribution},
                            marg_in1::Any,
                            marg_slope::Any,
                            marg_offset::Any,
                            marg_noise::GammaDistribution,
                            marg_out::Any)
    # Variational update function, takes the marginals as input.
    # Sends to any interface carrying a Gaussian message, while using precision parameterized noise
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Ensure right parameterization
    for param = [marg_in1, marg_slope, marg_offset, marg_out] # just get all of them, all marginals need to be defined anyway
        (param == nothing) || ensureMWParametrization!(param)
    end
    # Get the variables beforehand for more readable update equations
    if marg_slope != nothing
        mu_a = marg_slope.m[1]
        gam_a = marg_slope.W[1, 1] # gamma_a
    end
    if marg_offset != nothing
        mu_b = marg_offset.m[1]
        gam_b = marg_offset.W[1, 1] # gamma_b
    end
    if marg_out != nothing
        mu_x2 = marg_out.m[1]
        gam_x2 = marg_out.W[1, 1] # gamma_x2
    end
    if marg_in1 != nothing
        mu_x1 = marg_in1.m[1]
        gam_x1 = marg_in1.W[1, 1] # gamma_x1
    end
    a_gam = marg_noise.a
    b_gam = marg_noise.b

    if outbound_interface_id == 1 # in1
        dist_out.m = [(mu_a*(mu_x2 - mu_b))/(inv(gam_a) + mu_a^2)]
        dist_out.V = nothing
        dist_out.W = reshape([(a_gam*(mu_a^2 + inv(gam_a)))/(b_gam)], 1, 1)
        dist_out.xi = nothing
    elseif outbound_interface_id == 2 # a
        dist_out.m = [(mu_x1*(mu_x2 - mu_b))/(inv(gam_x1) + mu_x1^2)]
        dist_out.V = nothing
        dist_out.W = reshape([(a_gam*(mu_x1^2 + inv(gam_x1)))/b_gam], 1, 1)
        dist_out.xi = nothing
    elseif outbound_interface_id == 3 # b
        dist_out.m = [mu_x2 - mu_a*mu_x1]
        dist_out.V = nothing
        dist_out.W = reshape([a_gam/b_gam], 1, 1)
        dist_out.xi = nothing
    elseif outbound_interface_id == 5 # out
        dist_out.m = [mu_a*mu_x1 + mu_b]
        dist_out.V = nothing
        dist_out.W = reshape([a_gam/b_gam], 1, 1)
        dist_out.xi = nothing
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GaussianDistribution},
                            marg_in1::Any,
                            marg_slope::Any,
                            marg_offset::Any,
                            marg_noise::InverseGammaDistribution,
                            marg_out::Any)
    # Variational update function, takes the marginals as input.
    # Sends to any interface carrying a Gaussian message, while using variance parameterized noise
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Ensure right parameterization
    for param = [marg_in1, marg_slope, marg_offset, marg_out] # just get all of them, all marginals need to be defined anyway
        (param==nothing) || ensureMVParametrization!(param)
    end
    # Get the variables beforehand for more readable update equations
    if marg_slope != nothing
        mu_a = marg_slope.m[1]
        s_a = marg_slope.V[1, 1] # sigma_a^2
    end
    if marg_offset != nothing
        mu_b = marg_offset.m[1]
        s_b = marg_offset.V[1, 1] # sigma_b^2
    end
    if marg_out != nothing
        mu_x2 = marg_out.m[1]
        s_x2 = marg_out.V[1, 1] # sigma_x2^2
    end
    if marg_in1 != nothing
        mu_x1 = marg_in1.m[1]
        s_x1 = marg_in1.V[1, 1] # sigma_x1^2
    end
    a_s = marg_noise.a
    b_s = marg_noise.b

    if outbound_interface_id == 1 # in1
        dist_out.m = [(mu_a*(mu_x2 - mu_b))/(s_a + mu_a^2)]
        dist_out.V = reshape([(a_s + 1)/(b_s*(mu_a^2 + s_a))], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    elseif outbound_interface_id == 2 # a
        dist_out.m = [(mu_x1*(mu_x2 - mu_b))/(s_x1 + mu_x1^2)]
        dist_out.V = reshape([(a_s + 1)/(b_s*(mu_x1^2 + s_x1))], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    elseif outbound_interface_id == 3 # b
        dist_out.m = [mu_x2 - mu_a*mu_x1]
        dist_out.V = reshape([(a_s + 1)/b_s], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    elseif outbound_interface_id == 5 # out
        dist_out.m = [mu_a*mu_x1 + mu_b]
        dist_out.V = reshape([(a_s + 1)/b_s], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{GammaDistribution},
                            marg_in1::GaussianDistribution,
                            marg_slope::GaussianDistribution,
                            marg_offset::GaussianDistribution,
                            ::Any,
                            marg_out::GaussianDistribution)
    # Variational update function, takes the marginals as input.
    # Sends to precision parameterized noise.
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Ensure right parameterization
    for param = [marg_in1, marg_slope, marg_offset, marg_out] # just get all of them, all marginals need to be defined anyway
        ensureMWParametrization!(param)
    end

    # Get the variables beforehand for more readable update equations
    mu_a = marg_slope.m[1]
    gam_a = marg_slope.W[1, 1] # gamma_a
    mu_b = marg_offset.m[1]
    gam_b = marg_offset.W[1, 1] # gamma_b
    mu_x2 = marg_out.m[1]
    gam_x2 = marg_out.W[1, 1] # gamma_x2
    mu_x1 = marg_in1.m[1]
    gam_x1 = marg_in1.W[1, 1] # gamma_x1

    if outbound_interface_id == 4 # noise_N
        dist_out.a = 1.5
        dist_out.b = 0.5*((mu_x2 - mu_a*mu_x1 - mu_b)^2 - (mu_a^2)*(mu_x1^2) + (mu_a^2 + inv(gam_a))*(mu_x1^2 + inv(gam_x1)) + inv(gam_b) + inv(gam_x2))
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_value_type::Type{InverseGammaDistribution},
                            marg_in1::GaussianDistribution,
                            marg_slope::GaussianDistribution,
                            marg_offset::GaussianDistribution,
                            ::Any,
                            marg_out::GaussianDistribution)
    # Variational update function, takes the marginals as input.
    # Sends to variance parameterized noise.
    # Derivation for the update rule can be found in the derivations notebook.

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_value_type).value

    # Ensure right parameterization
    for param = [marg_in1, marg_slope, marg_offset, marg_out] # just get all of them, all marginals need to be defined anyway
        ensureMVParametrization!(param)
    end

    # Get the variables beforehand for more readable update equations
    mu_a = marg_slope.m[1]
    s_a = marg_slope.V[1, 1] # sigma_a^2
    mu_b = marg_offset.m[1]
    s_b = marg_offset.V[1, 1] # sigma_b^2
    mu_x2 = marg_out.m[1]
    s_x2 = marg_out.V[1, 1] # sigma_x2^2
    mu_x1 = marg_in1.m[1]
    s_x1 = marg_in1.V[1, 1] # sigma_x1^2

    if outbound_interface_id == 4 # s_N
        dist_out.a = -0.5
        dist_out.b = 0.5*((mu_x2 - mu_a*mu_x1 - mu_b)^2 - (mu_a^2)*(mu_x1^2) + (mu_a^2 + s_a)*(mu_x1^2 + s_x1) + s_b + s_x2)
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    # Set the outbound message
    return node.interfaces[outbound_interface_id].message
end