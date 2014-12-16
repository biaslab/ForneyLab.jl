############################################
# GainAdditionCompositeNode
############################################
# Description:
#   Composite gain-addition node.
#
#            | in1
#            |
#        ____|____
#        |   v   |
#        |  [A]  |
#        |   |   |
#    in2 |   v   | out
#   -----|->[+]--|---->
#        |_______|
#
#
#   out = A*in1 + in2
#   Example:
#       GainAdditionCompositeNode([1.0]; name="my_node")
#   Gain A is fixed and defined through the constructor.
#
# Interface ids, (names) and supported message types:
#   1. in1:
#       Message{GaussianDistribution}
#   2. in2:
#       Message{GaussianDistribution}
#   3. out:
#       Message{GaussianDistribution}
#
############################################
export GainAdditionCompositeNode

type GainAdditionCompositeNode <: CompositeNode
    # Basic node properties.
    # use_composite_update_rules is a flag for indicating use of
    # shortcut update rules. When set to true (default), the outbound
    # message is calculated via the provided update rules.
    # When set to false, messages are passed through the internal
    # graph of the composite node.
    A::Array
    use_composite_update_rules::Bool
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Pointers to internal nodes
    addition_node::AdditionNode
    fixed_gain_node::FixedGainNode
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface

    function GainAdditionCompositeNode(A::Union(Array{Float64},Float64)=1.0, use_composite_update_rules::Bool=true; name=unnamedStr())
        if typeof(A)==Float64
            A = fill!(Array(Float64,1,1),A)
        elseif use_composite_update_rules
            # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
            # In case we don't use composite update rules, A is passed to the internal FixedGainNode.
            A = ensureMatrix(deepcopy(A))
        end
        self = new(A, use_composite_update_rules, name, Array(Interface, 3))

        # Define the internals of the composite node
        self.addition_node = AdditionNode(name="$(name)_internal_addition")
        self.fixed_gain_node = FixedGainNode(A, name="$(name)_internal_gain")
        Edge(self.fixed_gain_node.out, self.addition_node.in1, GaussianDistribution, add_to_graph=false) # Internal edge

        named_handle_list = [:in1, :in2, :out]
        for i = 1:length(named_handle_list)
            self.interfaces[i] = Interface(self) # Initialize the composite node interfaces belonging to the composite node itself.
            setfield!(self, named_handle_list[i], self.interfaces[i]) # Init named interface handles
        end

        # Initialize the interfaces as references to the internal node interfaces.
        self.interfaces[1].child = self.fixed_gain_node.in1
        self.interfaces[2].child = self.addition_node.in2
        self.interfaces[3].child = self.addition_node.out

        # Set internal message passing schedules
        self.in1.internal_schedule = convert_to_schedule([self.addition_node.in1, self.fixed_gain_node.in1])
        self.in2.internal_schedule = convert_to_schedule([self.fixed_gain_node.out, self.addition_node.in2])
        self.out.internal_schedule = convert_to_schedule([self.fixed_gain_node.out, self.addition_node.out])

        return self
    end
end

isDeterministic(::GainAdditionCompositeNode) = true

############################################
# GaussianDistribution methods
############################################

# Rule set for forward propagation ({in1,in2}-->out)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainAdditionMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, m_y::Array{T, 1}) = m_x + A*m_y
forwardGainAdditionVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x + A*V_y*A'
forwardGainAdditionWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x - W_x * A * inv(W_y+A'*W_x*A) * A' * W_x
forwardGainAdditionXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = xi_x + W_x*A*inv(W_y+A'*W_x*A)*(xi_y-A'*xi_x)

# Rule set for backward propagation ({in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardIn2GainAdditionMRule{T<:Number}(A::Array{T, 2}, m_y::Array{T, 1}, m_z::Array{T, 1}) = m_z - A*m_y
backwardIn2GainAdditionVRule{T<:Number}(A::Array{T, 2}, V_y::Array{T, 2}, V_z::Array{T, 2}) = V_z + A*V_y*A'
backwardIn2GainAdditionWRule{T<:Number}(A::Array{T, 2}, W_y::Array{T, 2}, W_z::Array{T, 2}) = W_z - W_z * A * inv(W_y+A'*W_z*A) * A' * W_z
backwardIn2GainAdditionXiRule{T<:Number}(A::Array{T, 2}, xi_y::Array{T, 1}, xi_z::Array{T, 1}, W_y::Array{T, 2}, W_z::Array{T, 2}) = xi_z - W_z*A*inv(W_y+A'*W_z*A)*(xi_y+A'*xi_z)

# Forward to OUT
function updateNodeMessage!(node::GainAdditionCompositeNode,
                            outbound_interface_id::Int,
                            in1::Message{GaussianDistribution},
                            in2::Message{GaussianDistribution},
                            ::Nothing)

    if outbound_interface_id == 3
        if !node.use_composite_update_rules
            node.interfaces[outbound_interface_id].message = executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
        else
            dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload
    
            dist_1 = in1.payload
            dist_2 = in2.payload

            # Select parameterization
            # Order is from least to most computationally intensive
            if dist_1.m != nothing && dist_1.V != nothing && dist_2.m != nothing && dist_2.V != nothing
                dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
                dist_out.V  = forwardGainAdditionVRule(node.A, dist_2.V, dist_1.V)
                dist_out.W  = nothing
                dist_out.xi = nothing
            elseif dist_1.m != nothing && dist_1.W != nothing && dist_2.m != nothing && dist_2.W != nothing
                dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
                dist_out.V  = nothing
                dist_out.W  = forwardGainAdditionWRule(node.A, dist_2.W, dist_1.W)
                dist_out.xi = nothing
            elseif dist_1.xi != nothing && dist_1.W != nothing && dist_2.xi != nothing && dist_2.W != nothing
                dist_out.m  = nothing
                dist_out.V  = nothing
                dist_out.W  = forwardGainAdditionWRule(node.A, dist_2.W, dist_1.W)
                dist_out.xi = forwardGainAdditionXiRule(node.A, dist_2.xi, dist_1.xi, dist_2.W, dist_1.W)
            elseif (dist_1.m != nothing && dist_1.V != nothing) || (dist_2.m != nothing && dist_2.V != nothing)
                # Fallback: at least one inbound msg is in (m,V) parametrization
                # Convert the other one to (m,V)
                ensureMVParametrization!(dist_1)
                ensureMVParametrization!(dist_2)
                dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
                dist_out.V  = forwardGainAdditionVRule(node.A, dist_2.V, dist_1.V)
                dist_out.W  = nothing
                dist_out.xi = nothing
            elseif (dist_1.m != nothing && dist_1.W != nothing) || (dist_2.m != nothing && dist_2.W != nothing)
                # Fallback: at least one inbound msg is in (m,W) parametrization
                # Convert the other one to (m,W)
                ensureMWParametrization!(dist_1)
                ensureMWParametrization!(dist_2)
                dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
                dist_out.V  = nothing
                dist_out.W  = forwardGainAdditionWRule(node.A, dist_2.W, dist_1.W)
                dist_out.xi = nothing
            else
                # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
                ensureMVParametrization!(dist_1)
                ensureMVParametrization!(dist_2)
                dist_out.m  = forwardGainAdditionMRule(node.A, dist_2.m, dist_1.m)
                dist_out.V  = forwardGainAdditionVRule(node.A, dist_2.V, dist_1.V)
                dist_out.W  = nothing
                dist_out.xi = nothing
            end
        end
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    return node.interfaces[outbound_interface_id].message
end

# Backward to IN2
function updateNodeMessage!(node::GainAdditionCompositeNode,
                            outbound_interface_id::Int,
                            in1::Message{GaussianDistribution},
                            ::Nothing,
                            out::Message{GaussianDistribution})

    if outbound_interface_id == 2
        if !node.use_composite_update_rules
            node.interfaces[outbound_interface_id].message = executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
        else
            dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload

            dist_1 = in1.payload
            dist_3 = out.payload

            # Select parameterization
            # Order is from least to most computationally intensive
            if dist_1.m != nothing && dist_1.V != nothing && dist_3.m != nothing && dist_3.V != nothing
                dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
                dist_out.V  = backwardIn2GainAdditionVRule(node.A, dist_1.V, dist_3.V)
                dist_out.W  = nothing
                dist_out.xi = nothing
            elseif dist_1.m != nothing && dist_1.W != nothing && dist_3.m != nothing && dist_3.W != nothing
                dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
                dist_out.V  = nothing
                dist_out.W  = backwardIn2GainAdditionWRule(node.A, dist_1.W, dist_3.W)
                dist_out.xi = nothing
            elseif dist_1.xi != nothing && dist_1.W != nothing && dist_3.xi != nothing && dist_3.W != nothing
                dist_out.m  = nothing
                dist_out.V  = nothing
                dist_out.W  = backwardIn2GainAdditionWRule(node.A, dist_1.W, dist_3.W)
                dist_out.xi = backwardIn2GainAdditionXiRule(node.A, dist_1.xi, dist_3.xi, dist_1.W, dist_3.W)
            elseif (dist_1.m != nothing && dist_1.V != nothing) || (dist_3.m != nothing && dist_3.V != nothing)
                # Fallback: at least one inbound msg is in (m,V) parametrization
                # Convert the other one to (m,V)
                ensureMVParametrization!(dist_1)
                ensureMVParametrization!(dist_3)
                dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
                dist_out.V  = backwardIn2GainAdditionVRule(node.A, dist_1.V, dist_3.V)
                dist_out.W  = nothing
                dist_out.xi = nothing
            elseif (dist_1.m != nothing && dist_1.W != nothing) || (dist_3.m != nothing && dist_3.W != nothing)
                # Fallback: at least one inbound msg is in (m,W) parametrization
                # Convert the other one to (m,W)
                ensureMWParametrization!(dist_1)
                ensureMWParametrization!(dist_3)
                dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
                dist_out.V  = nothing
                dist_out.W  = backwardIn2GainAdditionWRule(node.A, dist_1.W, dist_3.W)
                dist_out.xi = nothing
            else
                # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
                ensureMVParametrization!(dist_1)
                ensureMVParametrization!(dist_3)
                dist_out.m  = backwardIn2GainAdditionMRule(node.A, dist_1.m, dist_3.m)
                dist_out.V  = backwardIn2GainAdditionVRule(node.A, dist_1.V, dist_3.V)
                dist_out.W  = nothing
                dist_out.xi = nothing
            end
        end
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    return node.interfaces[outbound_interface_id].message
end

# Backward to IN1
function updateNodeMessage!(node::GainAdditionCompositeNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            in2::Message{GaussianDistribution},
                            out::Message{GaussianDistribution})

    if outbound_interface_id == 1
        # We don't have a shortcut rule for this one, so we use the internal nodes to calculate the outbound msg
        node.interfaces[outbound_interface_id].message = executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

    return node.interfaces[outbound_interface_id].message
end