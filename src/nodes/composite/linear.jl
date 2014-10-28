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
#       Message{GaussianDistribution}
#       GaussianDistribution (marginal)
#   2. slope:
#       Message{Float64}
#       GaussianDistribution (marginal)
#   3. offset:
#       Message{Float64}
#       GaussianDistribution (marginal)
#   4. noise:
#       Message{Float64}
#       InvertedGammaDistribution (marginal)
#       GammaDistribution (marginal)
#   5. out:
#       Message{GaussianDistribution}
#       GaussianDistribution (marginal)
#
############################################

export LinearCompositeNode

type LinearCompositeNode <: CompositeNode
    # Basic node properties.
    # use_composite_update_rules is a flag for indicating use of shortcut update rules.
    # The linear composite node only has explicit variational update rules and thus does not need an explicit internal structure.

    use_composite_update_rules::Bool
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Field for form checking
    form::ASCIIString
    # Helper fields filled by constructor
    in1::Interface
    slope::Interface
    offset::Interface
    noise::Interface
    out::Interface

    function LinearCompositeNode(use_composite_update_rules::Bool=true; name = "unnamed", form::ASCIIString="moment")
        if use_composite_update_rules == false # Check
            error("LinearCompositeNode $(name) does not support explicit internal message passing")
        end

        self = new(use_composite_update_rules, name, Array(Interface, 5), form)

        # Set up the interfaces
        param_list = [:in1, :slope, :offset, :noise, :out]
        for i = 1:length(param_list)
            self.interfaces[i] = Interface(self) # Construct interface
            setfield!(self, param_list[i], self.interfaces[i]) # Set named interfaces
        end

        return self
    end
end


############################################
# Standard update functions
############################################

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            msg_in1::Message{GaussianDistribution},
                            msg_slope::Message{Float64},
                            msg_offset::Message{Float64},
                            msg_noise::Message{Float64},
                            ::Nothing)
    # Forward message to out, same rules as GainAdditionCompositeNode from Korl
    node.form == "precision" || error("You need to specify the 'precision' form when constructing $(typeof(node)) $(node.name) in order to work with mean-precision parametrization")

    ensureMWParametrization!(msg_in1.payload)

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload
    dist_out.m = forwardGainAdditionMRule(msg_slope.payload*eye(1), [msg_offset.payload], msg_in1.payload.m)
    dist_out.W = forwardGainAdditionWRule(msg_slope.payload*eye(1), msg_noise.payload*eye(1), msg_in1.payload.W)
    dist_out.V = nothing
    dist_out.xi = nothing

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            msg_slope::Message{Float64},
                            msg_offset::Message{Float64},
                            msg_noise::Message{Float64},
                            msg_out::Message{GaussianDistribution})
    # Backward message to in1, same rules as GainAdditionCompositeNode from Korl
    node.form == "precision" || error("You need to specify the 'precision' form when constructing $(typeof(node)) $(node.name) in order to work with mean-precision parametrization")

    ensureMWParametrization!(msg_out.payload)

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # No explicit backward rule exists, so stack them manually
    m_y = backwardAdditionMRule([msg_offset.payload], msg_out.payload.m) 
    W_y = backwardAdditionWRule(msg_noise.payload*eye(1), msg_out.payload.W)
    dist_out.m = backwardFixedGainMRule(inv(msg_slope.payload*eye(1)), m_y)
    dist_out.W = backwardFixedGainWRule(msg_slope.payload*eye(1), W_y)
    dist_out.V = nothing
    dist_out.xi = nothing

    return node.interfaces[outbound_interface_id].message
end


############################################
# Variational update functions
############################################

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            dist_gam::GammaDistribution,
                            dist_y::GaussianDistribution)
    # Backward message over x
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v   Q(gam)~Gam         
    #  ---->[  L  ]<----
    #        ^   |
    #      | |   | Q(y)~N
    #      v |   v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMWParametrization!(dist_a)
    ensureMDefined!(dist_b)
    ensureMDefined!(dist_y)
    
    mu_a = dist_a.m[1]
    gam_a = dist_a.W[1,1]
    mu_y = dist_y.m[1]
    mu_b = dist_b.m[1]
    a_gam = dist_gam.a
    b_gam = dist_gam.b

    if is(node.interfaces[outbound_interface_id], node.in1)
        dist_out.m = [(mu_a*(mu_y - mu_b))/(inv(gam_a) + mu_a^2)]
        dist_out.V = nothing
        dist_out.W = reshape([(a_gam*(mu_a^2 + inv(gam_a)))/(b_gam)], 1, 1)
        dist_out.xi = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            dist_x::GaussianDistribution,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            dist_gam::GammaDistribution,
                            ::Nothing)
    # Forward message over y
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v   Q(gam)~Gam         
    #  ---->[  L  ]<----
    #        ^   |
    # Q(x)~N |   | |
    #        |   v v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMDefined!(dist_a)
    ensureMDefined!(dist_b)
    ensureMDefined!(dist_x)

    mu_a = dist_a.m[1]
    mu_b = dist_b.m[1]
    mu_x = dist_x.m[1]
    a_gam = dist_gam.a
    b_gam = dist_gam.b

    if is(node.interfaces[outbound_interface_id], node.out)
        dist_out.m = [mu_a*mu_x + mu_b]
        dist_out.V = nothing
        dist_out.W = reshape([a_gam/b_gam], 1, 1)
        dist_out.xi = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end


function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            dist_s::InverseGammaDistribution,
                            dist_y::GaussianDistribution)
    # Backward message over x
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v   Q(s)~Ig         
    #  ---->[  L  ]<----
    #        ^   |
    #      | |   | Q(y)~N
    #      v |   v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMVParametrization!(dist_a)
    ensureMDefined!(dist_b)
    ensureMDefined!(dist_y)
    
    mu_a = dist_a.m[1]
    s_a = dist_a.V[1,1]
    mu_y = dist_y.m[1]
    mu_b = dist_b.m[1]
    a_s = dist_s.a
    b_s = dist_s.b

    if is(node.interfaces[outbound_interface_id], node.in1)
        dist_out.m = [(mu_a*(mu_y - mu_b))/(s_a + mu_a^2)]
        dist_out.V = reshape([(a_s + 1)/(b_s*(mu_a^2 + s_a))], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            dist_x::GaussianDistribution,
                            ::Nothing,
                            dist_b::GaussianDistribution,
                            dist_s::InverseGammaDistribution,
                            dist_y::GaussianDistribution)
    # Backward message over slope
    #
    #        Q(b)~N
    #          |
    #   <--    v   Q(s)~Ig         
    #  ---->[  L  ]<----
    #        ^   |
    #        |   v
    #   Q(x)~N   Q(y)~N

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMVParametrization!(dist_x)
    ensureMDefined!(dist_y)
    ensureMDefined!(dist_b)

    mu_x = dist_x.m[1]
    s_x = dist_x.V[1,1]
    mu_y = dist_y.m[1]
    mu_b = dist_b.m[1]
    a_s = dist_s.a
    b_s = dist_s.b

    if is(node.interfaces[outbound_interface_id], node.slope)
        dist_out.m = [(mu_x*(mu_y - mu_b))/(s_x + mu_x^2)]
        dist_out.V = reshape([(a_s + 1)/(b_s*(mu_x^2 + s_x))], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            dist_x::GaussianDistribution,
                            dist_a::GaussianDistribution,
                            ::Nothing,
                            dist_s::InverseGammaDistribution,
                            dist_y::GaussianDistribution)
    # Backward message over offset
    #
    #        ^ |
    #        | |
    # Q(b)~N   v   Q(s)~Ig         
    #  ---->[  L  ]<----
    #        ^   |
    #        |   v
    #   Q(x)~N   Q(y)~N

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMDefined!(dist_x)
    ensureMDefined!(dist_y)
    ensureMDefined!(dist_a)

    mu_x = dist_x.m[1]
    mu_a = dist_a.m[1]
    mu_y = dist_y.m[1]
    a_s = dist_s.a
    b_s = dist_s.b

    if is(node.interfaces[outbound_interface_id], node.offset)
        dist_out.m = [mu_y - mu_a*mu_x]
        dist_out.V = reshape([(a_s + 1)/b_s], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{InverseGammaDistribution},
                            dist_x::GaussianDistribution,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            ::Nothing,
                            dist_y::GaussianDistribution)
    # Backward message over noise
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v    -->         
    #  ---->[  L  ]<----
    #        ^   |
    #        |   v
    #   Q(x)~N   Q(y)~N

    (node.form == "moment") || error("LinearCompositeNode $(node.name) must be declared in moment form")
    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMVParametrization!(dist_x)
    ensureMVParametrization!(dist_a)
    ensureMVParametrization!(dist_b)
    ensureMVParametrization!(dist_y)

    mu_a = dist_a.m[1]
    s_a = dist_a.V[1,1]
    mu_b = dist_b.m[1]
    s_b = dist_b.V[1,1]
    mu_y = dist_y.m[1]
    s_y = dist_y.V[1,1]
    mu_x = dist_x.m[1]
    s_x = dist_x.V[1,1]

    if is(node.interfaces[outbound_interface_id], node.noise)
        dist_out.a = -0.5
        dist_out.b = 0.5*((mu_y - mu_a*mu_x - mu_b)^2 - (mu_a^2)*(mu_x^2) + (mu_a^2 + s_a)*(mu_x^2 + s_x) + s_b + s_y)
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            dist_x::GaussianDistribution,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            dist_s::InverseGammaDistribution,
                            ::Nothing)
    # Forward message over y
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v   Q(s)~Ig         
    #  ---->[  L  ]<----
    #        ^   |
    # Q(x)~N |   | |
    #        |   v v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMDefined!(dist_a)
    ensureMDefined!(dist_b)
    ensureMDefined!(dist_x)

    mu_a = dist_a.m[1]
    mu_b = dist_b.m[1]
    mu_x = dist_x.m[1]
    a_s = dist_s.a
    b_s = dist_s.b

    if is(node.interfaces[outbound_interface_id], node.out)
        dist_out.m = [mu_a*mu_x + mu_b]
        dist_out.V = reshape([(a_s + 1)/b_s], 1, 1)
        dist_out.W = nothing
        dist_out.xi = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end


############################################
# (Structured) variational update functions
############################################

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Nothing,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            dist_gam::GammaDistribution,
                            ::Message{GaussianDistribution})
    # Backward message over x
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v   Q(gam)~Gam         
    #  ---->[  L  ]<----
    #        ^   |
    #      | |   | y~N
    #      v |   v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMVParametrization!(dist_a)
    ensureMDefined!(dist_b)

    mu_a = dist_a.m[1]
    sig_a_2 = dist_a.V[1,1]
    mu_b = dist_b.m[1]
    gam_a = dist_gam.a
    gam_b = dist_gam.b

    if is(node.interfaces[outbound_interface_id], node.in1)
        dist_out.m = [(mu_a*mu_b)/(sig_a_2+mu_a^2)]
        dist_out.W = reshape([(gam_a*(sig_a_2 + mu_a^2))/(gam_b)],1,1)
        dist_out.xi = nothing
        dist_out.V = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end   

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            dist_x::GaussianDistribution,
                            ::Nothing,
                            dist_b::GaussianDistribution,
                            dist_gam::GammaDistribution,
                            dist_y::GaussianDistribution)
    if is(dist_x, dist_y)
        # Structured

        # Backward message over slope
        #
        #        Q(b)~N
        #          |
        #   <--    v   Q(gam)~Gam         
        #  ---->[  L  ]<----
        #        ^   |
        #        |   v
        #       Q(x,y)~N

        dist_xy = dist_x
        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

        # Ensure right parameterization
        ensureMVParametrization!(dist_xy)
        ensureMDefined!(dist_b)

        a_gam = dist_gam.a
        b_gam = dist_gam.b
        mu_b = dist_b.m[1]
        mu_x = dist_xy.m[1]
        mu_y = dist_xy.m[2]
        rho_xx = dist_xy.V[1,1]
        rho_xy = dist_xy.V[1,2]

        if is(node.interfaces[outbound_interface_id], node.slope)
            dist_out.m = [(mu_x*mu_y + rho_xy - mu_b*mu_x)/(mu_x^2 + rho_xx)]
            dist_out.W = reshape([(a_gam*(mu_x^2+rho_xx))/(b_gam)],1,1)
            dist_out.xi = nothing
            dist_out.V = nothing
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    else
        # Naive

        # Backward message over slope
        #
        #        Q(b)~N
        #          |
        #   <--    v   Q(gam)~Gam         
        #  ---->[  L  ]<----
        #        ^   |
        #        |   v
        #   Q(x)~N   Q(y)~N

        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

        # Ensure right parameterization
        ensureMWParametrization!(dist_x)
        ensureMDefined!(dist_y)
        ensureMDefined!(dist_b)

        mu_x = dist_x.m[1]
        gam_x = dist_x.W[1,1]
        mu_y = dist_y.m[1]
        mu_b = dist_b.m[1]
        a_gam = dist_gam.a
        b_gam = dist_gam.b

        if is(node.interfaces[outbound_interface_id], node.slope)
            dist_out.m = [(mu_x*(mu_y - mu_b))/(inv(gam_x) + mu_x^2)]
            dist_out.V = nothing
            dist_out.W = reshape([(a_gam*(mu_x^2 + inv(gam_x)))/b_gam], 1, 1)
            dist_out.xi = nothing
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            dist_x::GaussianDistribution,
                            dist_a::GaussianDistribution,
                            ::Nothing,
                            dist_gam::GammaDistribution,
                            dist_y::GaussianDistribution)
    if is(dist_x, dist_y)
        # Structured

        # Backward message over slope
        #
        #        ^ |
        #        | |
        # Q(a)~N   v   Q(gam)~Gam         
        #  ---->[  L  ]<----
        #        ^   |
        #        |   v
        #       Q(x,y)~N

        dist_xy = dist_x
        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

        # Ensure right parameterization
        ensureMVParametrization!(dist_xy)
        ensureMDefined!(dist_a)

        a_gam = dist_gam.a
        b_gam = dist_gam.b
        mu_a = dist_a.m[1]
        mu_x = dist_xy.m[1]
        mu_y = dist_xy.m[2]

        if is(node.interfaces[outbound_interface_id], node.offset)
            dist_out.m = [mu_y - mu_a*mu_x]
            dist_out.W = reshape([a_gam/b_gam],1,1)
            dist_out.xi = nothing
            dist_out.V = nothing
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    else
        # Naive

        # Backward message over offset
        #
        #        ^ |
        #        | |
        # Q(b)~N   v   Q(gam)~Gam         
        #  ---->[  L  ]<----
        #        ^   |
        #        |   v
        #   Q(x)~N   Q(y)~N

        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

        # Ensure right parameterization
        ensureMDefined!(dist_x)
        ensureMDefined!(dist_y)
        ensureMDefined!(dist_a)

        mu_x = dist_x.m[1]
        mu_a = dist_a.m[1]
        mu_y = dist_y.m[1]
        a_gam = dist_gam.a
        b_gam = dist_gam.b

        if is(node.interfaces[outbound_interface_id], node.offset)
            dist_out.m = [mu_y - mu_a*mu_x]
            dist_out.V = nothing
            dist_out.W = reshape([a_gam/b_gam], 1, 1)
            dist_out.xi = nothing
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GammaDistribution},
                            dist_x::GaussianDistribution,
                            dist_a::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            ::Nothing,
                            dist_y::GaussianDistribution)
    if is(dist_x, dist_y)
        # Structured

        # Backward message over noise
        #
        #        Q(b)~N
        #          |
        # Q(a)~N   v    -->         
        #  ---->[  L  ]<----
        #        ^   |
        #        |   v
        #       Q(x,y)~N

        dist_xy = dist_x
        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

        # Ensure right parameterization
        ensureMVParametrization!(dist_xy)
        ensureMVParametrization!(dist_a)
        ensureMVParametrization!(dist_b)

        mu_a = dist_a.m[1]
        sig_a_2 = dist_a.V[1,1]
        mu_b = dist_b.m[1]
        sig_b_2 = dist_b.V[1,1]
        mu_x = dist_xy.m[1]
        mu_y = dist_xy.m[2]
        rho_xx = dist_xy.V[1,1]
        rho_xy = dist_xy.V[1,2]
        rho_yy = dist_xy.V[2,2]

        if is(node.interfaces[outbound_interface_id], node.noise)
            dist_out.a = 1.5
            dist_out.b = 0.5*(mu_y^2*rho_yy - mu_a*(mu_x*mu_y + rho_xy) - mu_b*mu_y + (mu_a^2 + sig_a_2)*(mu_x^2 + rho_xx) + mu_a*mu_b*mu_x + mu_b^2 + sig_b_2)
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    else
        # Naive

        # Backward message over noise
        #
        #        Q(b)~N
        #          |
        # Q(a)~N   v    -->         
        #  ---->[  L  ]<----
        #        ^   |
        #        |   v
        #   Q(x)~N   Q(y)~N

        (node.form == "precision") || error("LinearCompositeNode $(node.name) must be declared in precision form")
        dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

        # Ensure right parameterization
        ensureMWParametrization!(dist_x)
        ensureMWParametrization!(dist_a)
        ensureMWParametrization!(dist_b)
        ensureMWParametrization!(dist_y)

        mu_a = dist_a.m[1]
        gam_a = dist_a.W[1,1]
        mu_b = dist_b.m[1]
        gam_b = dist_b.W[1,1]
        mu_y = dist_y.m[1]
        gam_y = dist_y.W[1,1]
        mu_x = dist_x.m[1]
        gam_x = dist_x.W[1,1]

        if is(node.interfaces[outbound_interface_id], node.noise)
            dist_out.a = 1.5
            dist_out.b = 0.5*((mu_y - mu_a*mu_x - mu_b)^2 - (mu_a^2)*(mu_x^2) + (mu_a^2 + inv(gam_a))*(mu_x^2 + inv(gam_x)) + inv(gam_b) + inv(gam_y))
        else
            error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
        end
    end

    return node.interfaces[outbound_interface_id].message
end

function updateNodeMessage!(node::LinearCompositeNode,
                            outbound_interface_id::Int,
                            outbound_message_payload_type::Type{GaussianDistribution},
                            ::Message{GaussianDistribution},
                            ::GaussianDistribution,
                            dist_b::GaussianDistribution,
                            dist_gam::GammaDistribution,
                            ::Nothing)
    # Forward message over y
    #
    #        Q(b)~N
    #          |
    # Q(a)~N   v   Q(gam)~Gam         
    #  ---->[  L  ]<----
    #        ^   |
    #    x~N |   | |
    #        |   v v

    dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], outbound_message_payload_type).payload

    # Ensure right parameterization
    ensureMDefined!(dist_b)

    mu_b = dist_b.m[1]
    gam_a = dist_gam.a
    gam_b = dist_gam.b

     if is(node.interfaces[outbound_interface_id], node.out)
        dist_out.m = [mu_b]
        dist_out.W = reshape([gam_a/gam_b],1,1)
        dist_out.xi = nothing
        dist_out.V = nothing
    else
        error("Undefined inbound-outbound message type combination for node $(node.name) of type $(typeof(node)).")
    end

    return node.interfaces[outbound_interface_id].message
end 