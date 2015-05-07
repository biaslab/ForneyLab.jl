**************
 Nodes
**************

This chapter describes :ref:`node-anatomy`, :ref:`msg-calc-rules`, and :ref:`built-in-nodes`. For more details, have a look at the source files of the built-in nodes, such as ``src/nodes/addition.jl``. 

.. _node-anatomy:

The anatomy of nodes
--------------------

A factor graph node is always a subtype of ``abstract Node``. A node should contain at least the fields ``name`` and ``interfaces``. Let's look at the definition of the built-in :class:`AdditionNode`::

    type AdditionNode <: Node
        name::ASCIIString
        interfaces::Array{Interface,1}
        in1::Interface
        in2::Interface
        out::Interface

        function AdditionNode(; name=unnamedStr())
            self = new(name, Array(Interface, 3))

            named_handle_list = [:in1, :in2, :out]
            for i = 1:length(named_handle_list)
                self.interfaces[i] = Interface(self)
                setfield!(self, named_handle_list[i], self.interfaces[i])
            end

            return self
        end
    end

The fields ``in1``, ``in2``, and ``out`` are optional 'named handles' to make accessing the interfaces convenient. The ``interfaces`` array always contains one or more :class:`Interface` instances. The calling signature of a node constructor varies, but it always includes the optional keyword argument ``name``. 

.. type:: Interface
    
    An Interface belongs to a node and can be partnered to another ``Interface`` to form an :class:`Edge`. It can be viewed as an half-edge that can be combined with another half-edge to form a complete :class:`Edge`.
    ::

        type Interface
            node::Node # Reference to the Node it is part of
            edge::Union(AbstractEdge, Nothing) # Reference to the Edge it is part of
            partner::Union(Interface, Nothing) # Partner it is connected to
            child::Union(Interface, Nothing)   # For composite nodes
            message::Union(Message, Nothing)   # Outbound message
            internal_schedule::Array{Any, 1}   # For composite nodes: schedule to calc outbound msg
        end

    The ``message`` field may contain a :class:`Message`, which is the *outbound message* calculated according to the node function. This means that if an interface is the tail of an :class:`Edge`, its ``message`` field contains the *forward message* on that edge. Similarly, if the interface is the head of the edge, its ``message`` field contains the *backward message*. 

.. _msg-calc-rules:

Message calculation rules
-------------------------

A factor node captures a specific *node function*, which involves the variables that are represented by the connected edges. The :class:`AdditionNode` for example captures the addition function ``out = in1 + in2``. When running a message passing algorithm on a :class:`FactorGraph`, the node function specifies how the outbound messages are calculated from the inbound messages. An outbound message is calculated according to a *message calculation rule*. Message calculation rules are implemented for specific nodes and message types using multiple dispatch. 

ForneyLab supports the following message calculation rules:

.. function:: sumProduct!(node::Node, outbound_interface_id::Int, inbound_messages...)
    
    Calculates the outbound message from the incoming messages on the other interfaces according to the sum-product algorithm.
    Example implementation::

        function sumProduct!(node::AdditionNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Nothing)
            # Calculate outbound message on interface 3 (named "out") 
            # according to the sum-product rule

            # Save the calculated message on the interface

            # Return tuple ([calculation rule name]::Symbol, outbound_message::Message)
            return (:addition_gaussian_forward,
                    node.interfaces[outbound_interface_id].message)
        end

    The calling signature consists of:

    1. The node;
    2. The id (index in node.interfaces) of the outbound interface;
    3. The inbound messages on *all* interfaces of the node (ordered by interface id). The inbound message on the outbound inferface is always ``nothing``.

.. function:: vmp!(node::Node, outbound_interface_id::Int, marginals_and_messages...)

    Similar to :func:`sumProduct!`, but on some interfaces the approximate marginals are used instead of the incoming messages. This calculation rule is used for variational message passing (vmp).
    Example implementation::

        function vmp!(node::GaussianNode,
                            outbound_interface_id::Int,
                            ::Nothing,
                            marg_prec::GammaDistribution,
                            marg_y::GaussianDistribution)
            # Calculate outbound message on interface 1 (named "mean") 
            # according to the variational message passing rule

            # Save the calculated message on the interface

            # Return tuple ([calculation rule name]::Symbol, outbound_message::Message)
            return (:gaussian_backward_mean_gaussian_inverse_gamma,
                    node.interfaces[outbound_interface_id].message)
        end                            

    The calling signature consists of:

    1. The node;
    2. The id (index in node.interfaces) of the outbound interface;
    3. The messages or marginals on *all* interfaces of the node (ordered by interface id). The inbound message/marginal on the outbound inferface is always ``nothing``.

Not all message calculation rules have to be implemented for every node, just the ones that will be used. Similarly, the message calculation rule does not have to be implemented for a specific outbound interface of a specific node if that outbound message never has to be calculated.

.. _built-in-nodes:

Built-in nodes
--------------

The following built-in 'elementary' nodes are available in ForneyLab: :class:`AdditionNode`, :class:`EqualityNode`, :class:`ExponentialNode`, :class:`FixedGainNode`, :class:`GaussianNode`, :class:`TerminalNode`.

There are aso some built-in composite nodes: :class:`GainAdditionCompositeNode`, :class:`GainEqualityCompositeNode`.

