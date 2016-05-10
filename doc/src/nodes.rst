**************
 Nodes
**************

This chapter describes the implementation of nodes in ForneyLab. For more details on specific node types, use Julia's help functionality (``?AdditionNode``), or have a look at the source files, such as ``src/nodes/addition.jl``.


The anatomy of nodes
--------------------

Each node type is a subtype of abstract ``Node``. A node should contain at least the fields ``id``, ``interfaces`` and ``i``. Let's look at the definition of the built-in :class:`AdditionNode`::

    type AdditionNode <: Node
        id::Symbol
        interfaces::Array{Interface,1}
        i::Dict{Symbol, Interface}
    end

The field ``i`` stores 'named handles' to make accessing the interfaces convenient, for example if we want to access the output interface we can use ``node.i[:out]``. The ``interfaces`` array always contains one or more :class:`Interface` instances, and the index in this array is called the "interface id". The calling signature of a node constructor varies, but it always includes the optional keyword argument ``id``. A ``Node`` can be copied using ``copy(src::Node; id=:new_id)``, where ``:new_id`` will become the id of the copy. The copy contains the exact internal state of the original, but has no edges connected to it.

.. type:: Interface

    An Interface belongs to a node and can be partnered to another ``Interface`` to form an :class:`Edge`. It can be viewed as an half-edge that can be combined with another half-edge to form a complete :class:`Edge`.
    ::

        type Interface
            node::Node # Reference to the Node it is part of
            edge::Union(AbstractEdge, Void) # Reference to the Edge it is part of
            partner::Union(Interface, Void) # Partner it is connected to
            message::Union(Message, Void)   # Outbound message
        end

    The ``message`` field may contain a :class:`Message`, which is the *outbound message* calculated according to the node function. This means that if an interface is the tail of an :class:`Edge`, its ``message`` field contains the *forward message* on that edge. Similarly, if the interface is the head of the edge, its ``message`` field contains the *backward message*.


Message calculation rules
-------------------------

A factor node captures a specific *node function*, which involves the variables that are represented by the connected edges. The :class:`AdditionNode` for example captures the addition function: ``f(in1,in2,out) = δ(in1+in2-out)``. When running a message passing algorithm on a :class:`FactorGraph`, the node function specifies how the outbound messages are calculated from the inbound messages. An outbound message is calculated according to a message calculation rule, also simply called *rule*. Rules are implemented for specific nodes and message types, making heavy use of Julia's multiple dispatch system.

ForneyLab comes with the following built-in rules:

.. function:: sumProductRule!(node::Node, outbound_interface_index::Type{Val{i}}, outbound_dist::ProbabilityDistribution, inbound_messages...)

    Calculates the outbound message from the incoming messages on the other interfaces according to the sum-product algorithm.
    Example implementation::

        function sumProductRule!(node::AdditionNode,
                                 outbound_interface_index::Type{Val{3}},
                                 outbound_dist::Gaussian,
                                 msg_in1::Message{Gaussian},
                                 msg_in2::Message{Gaussian},
                                 msg_out::Any)
            # Calculate outbound message on interface 3 (:out)
            # according to the sum-product rule

            # Perform computations as in in-place operation on outbound_dist

            return outbound_dist
        end

    The calling signature consists of:

    1. The node;
    2. The interface id (index in node.interfaces) of the outbound interface as a value type;
    3. The payload of the message currently present on the outbound interface;
    4. The inbound messages on *all* interfaces of the node (ordered by interface id).

.. function:: variationalRule!(node::Node, outbound_interface_index::Type{Val{i}}, outbound_dist::ProbabilityDistribution, marginals_and_messages...)

    Similar to :func:`sumProductRule!`, but on some interfaces marginals are required instead of messages. This rule is used for variational message passing (vmp).
    Example implementation::

        function variationalRule!(node::GaussianNode,
                                  outbound_interface_index::Type{Val{1}},
                                  outbound_dist::Gaussian,
                                  marg_mean::Any,
                                  marg_prec::Gamma,
                                  marg_y::Gaussian)
            # Calculate outbound message on interface 1 (:mean)
            # according to the variational rule

            # Perform computations as in in-place operation on outbound_dist

            return outbound_dist
        end

    The calling signature consists of:

    1. The node;
    2. The interface id (index in node.interfaces) of the outbound interface as a value type;
    3. The payload of the message currently present on the outbound interface;
    4. The inbound messages/marginals on *all* interfaces of the node (ordered by interface id).

.. function:: expectationRule!(node::Node, outbound_interface_index::Type{Val{i}}, outbound_dist::Gaussian, inbound_messages...)

    Similar to :func:`sumProductRule!`, but also the inbound message on the outbound interface is consumed (this messages carries the cavity distrubution). This calculation rule is used in the expectation propagation algorithm.

    The calling signature consists of:

    1. The node;
    2. The interface id (index in node.interfaces) of the outbound interface as a value type;
    3. The payload of the message currently present on the outbound interface;
    4. The inbound messages on *all* interfaces of the node (ordered by interface id).

Not all message calculation rules have to be implemented for every node, just the ones that will be used. Similarly, the message calculation rule does not have to be implemented for a specific outbound interface of a specific node if that outbound message never has to be calculated.

To find out which message calculation rules are implemented for a specific node, use the ``rules`` function:

.. function:: rules(node_type::DataType, [rule::Function; outbound::Int=0])

    Print all message calculation rules implemented for ``node_type <: Node``.
    Optionally, the list can be restricted to a specific rule, such as ``sumProductRule!`` or ``variationalRule!``.
    If keyword argument ``outbound`` is passed, only the rules for that outbound interface id are listed.


Approximate rules
-----------------

A rule might not result in a convenient or tractable exact outbound message. In such cases, one might want to implement an *approximate rule*. For example, a difficult, non-Gaussian, outbound message might be approximated by a Gaussian outbound message using Laplace's method or by moment matching. An approximate rule should have an extra argument to specify the type of approximation. For example::

    function sumProductRule!(node::GaussianNode{:mean}{:precision},
                             outbound_interface_index::Type{Val{3}},
                             outbound_dist::Gaussian,
                             msg_mean::Message{Gaussian},
                             msg_prec::Message{Gamma},
                             msg_out::Any,
                             approx::Type{MomentMatching})
        # Approximate exact student's t message by a Gaussian through moment matching
        ...

        return outbound_dist
    end

The extra argument should be a subtype of ``ApproximationType``. Built-in approximation types are ``Laplace`` and ``MomentMatching``. By default, ForneyLab will only resort to an approximate rule if there is no exact rule that can handle the incoming messages. However, a the user can force a message to have a specific distribution type by passing the ``message_types`` keyword argument to the algorithm constructor. If there are multiple approximation types for the same approximating distribution type, the user can even specify the desired approximation::

    # Force the use of the Laplace approximation on my_interface
    msg_types = Dict{Interface,DataType}(
                    my_interface => Approximation{Gaussian, Laplace}
                )
    algo = SumProduct(message_types=msg_types)


Composite nodes
---------------

A node in which multiple node functions are combined into one node function is called a ``CompositeNode`` In ForneyLab 0.4 the composite node functionality is removed. It will be reintroduced in ForneyLab 0.5.
