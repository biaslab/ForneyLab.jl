**************
 Nodes
**************

This chapter describes :ref:`node-anatomy`, :ref:`msg-calc-rules`, :ref:`built-in-nodes`, and :ref:`composite-nodes`. For more details, have a look at the source files of the built-in nodes, such as ``src/nodes/addition.jl``.

.. _node-anatomy:

The anatomy of nodes
--------------------

A factor graph node is always a subtype of ``abstract Node``, and its name ends in "Node". A node should contain at least the fields ``name``, ``interfaces`` and ``i``. Let's look at the definition of the built-in :class:`AdditionNode`::

    type AdditionNode <: Node
        id::Symbol
        interfaces::Array{Interface,1}
        i::Dict{Symbol, Interface}
    end

The field ``i`` stores 'named handles' to make accessing the interfaces convenient, for example if we want to access the out interface we type ``node.i[:out]``. The ``interfaces`` array always contains one or more :class:`Interface` instances. The calling signature of a node constructor varies, but it always includes the optional keyword argument ``id``. A ``Node`` can be copied using ``copy(src::Node; id=:new_id)``, where ``:new_id`` will become the id of the copy. The copy contains the exact internal state of the original, but has no edges connected to it.

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

.. _msg-calc-rules:

Message calculation rules
-------------------------

A factor node captures a specific *node function*, which involves the variables that are represented by the connected edges. The :class:`AdditionNode` for example captures the addition function: ``f(in1,in2,out) = δ(out-in1-in2)``. When running a message passing algorithm on a :class:`FactorGraph`, the node function specifies how the outbound messages are calculated from the inbound messages. An outbound message is calculated according to a *message calculation rule*. Message calculation rules are implemented for specific nodes and message types using multiple dispatch.

ForneyLab supports the following message calculation rules:

.. function:: sumProduct!(node::Node, outbound_interface_index::Int, inbound_messages...)

    Calculates the outbound message from the incoming messages on the other interfaces according to the sum-product algorithm.
    Example implementation::

        function sumProduct!(node::AdditionNode,
                            outbound_interface_index::Int,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Void)
            # Calculate outbound message on interface 3 (:out)
            # according to the sum-product rule

            # Save the calculated message on the interface

            # Return tuple ([calculation rule name]::Symbol, outbound_message::Message)
            return (:addition_gaussian_forward,
                    node.interfaces[outbound_interface_index].message)
        end

    The calling signature consists of:

    1. The node;
    2. The index (index in node.interfaces) of the outbound interface;
    3. The inbound messages on *all* interfaces of the node (ordered by interface id). The inbound message on the outbound inferface is always ``nothing``.

.. function:: vmp!(node::Node, outbound_interface_index::Int, marginals_and_messages...)

    Similar to :func:`sumProduct!`, but on some interfaces the approximate marginals are used instead of the incoming messages. This calculation rule is used for variational message passing (vmp).
    Example implementation::

        function vmp!(node::GaussianNode,
                            outbound_interface_index::Int,
                            ::Void,
                            marg_prec::GammaDistribution,
                            marg_y::GaussianDistribution)
            # Calculate outbound message on interface 1 (:mean)
            # according to the variational message passing rule

            # Save the calculated message on the interface

            # Return tuple ([calculation rule name]::Symbol, outbound_message::Message)
            return (:gaussian_backward_mean_gaussian_inverse_gamma,
                    node.interfaces[outbound_interface_index].message)
        end

    The calling signature consists of:

    1. The node;
    2. The id (index in node.interfaces) of the outbound interface;
    3. The messages or marginals on *all* interfaces of the node (ordered by interface id). The inbound message/marginal on the outbound inferface is always ``nothing``.

.. function:: ep!(node::Node, outbound_interface_index::Int, inbound_messages...)

    Similar to :func:`sumProduct!`, but also the inbound message on the outbound interface is consumed (this messages carries the cavity distrubution). This calculation rule is used in the expectation propagation algorithm.

    The calling signature consists of:

    1. The node;
    2. The id (index in node.interfaces) of the outbound interface;
    3. The inbound messages on *all* interfaces of the node (ordered by interface id).

Not all message calculation rules have to be implemented for every node, just the ones that will be used. Similarly, the message calculation rule does not have to be implemented for a specific outbound interface of a specific node if that outbound message never has to be calculated.

.. _built-in-nodes:

Built-in nodes
--------------

The following built-in 'elementary' nodes are available in ForneyLab: :class:`AdditionNode`, :class:`EqualityNode`, :class:`ExponentialNode`, :class:`FixedGainNode`, :class:`GaussianNode`, :class:`SigmoidNode`, :class:`TerminalNode`.

There are also some built-in *combined nodes*, which combine two or more node functions into one for higher computational efficiency: :class:`GainAdditionNode`, :class:`GainEqualityNode`.

Elementary nodes
~~~~~~~~~~~~~~~~

.. type:: AdditionNode

    ::

               in2
               |
         in1   v  out
        ----->[+]----->

    :Node function: ``f(in1,in2,out) = δ(out-in1-in2)``
    :Interfaces:    1. ``i[:in1]``, 2. ``i[:in2]``, 3. ``i[:out]``
    :Construction:  ``AdditionNode(id="something")``

    Message computation rules:

    +-----------------+-----------------------------------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface                                      |
    + Rule            +-------------------------+-------------------------+-------------------------+
    |                 | 1                       | 2                       |  3                      |
    +=================+=========================+=========================+=========================+
    | sumProduct!     | ↓↑ ``Msg{Delta}``       | ↓↑ ``Msg{Delta}``       | ↓↑ ``Msg{Delta}``       |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↓↑ ``Msg{Gaussian}``    | ↓↑ ``Msg{Gaussian}``    | ↓↑ ``Msg{Gaussian}``    |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↑  ``Msg{Gaussian}``    | ↓  ``Msg{Gaussian}``    | ↓  ``Msg{Delta}``       |
    +                 +                         +-------------------------+-------------------------+
    |                 |                         | ↓  ``Msg{Delta}``       | ↓  ``Msg{Gaussian}``    |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↓  ``Msg{Gaussian}``    | ↑  ``Msg{Gaussian}``    | ↓  ``Msg{Delta}``       |
    +                 +-------------------------+                         +-------------------------+
    |                 | ↓  ``Msg{Delta}``       |                         | ↓  ``Msg{Gaussian}``    |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↓  ``Msg{Gaussian}``    | ↓  ``Msg{Gaussian}``    | ↑  ``Msg{Gaussian}``    |
    +                 +-------------------------+-------------------------+                         +
    |                 | ↓  ``Msg{Delta}``       | ↓  ``Msg{Delta}``       |                         |
    +-----------------+-------------------------+-------------------------+-------------------------+

.. type:: EqualityNode

    ::

               Y
               |
           X   v  Z
        ----->[=]----->

    :Node function: ``f(X,Y,Z) = δ(X-Z)δ(Y-Z)``
    :Interfaces:    1. ``i[1]``, 2. ``i[2]``, 3. ``i[3]``
    :Construction:  ``EqualityNode(id="something")``

    Message computation rules (\* = approximation):

    +-----------------+-----------------------------------------------------------------------------+
    |                 | Input/output (node is symmetrical in all interfaces)                        |
    + Rule            +-------------------------+---------------------------------------------------+
    |                 | Outbound interface      | Inbound interfaces                                |
    +=================+=========================+===================================================+
    | sumProduct!     | ``Msg{Delta}``          | ``Msg{Delta}`` and ``Msg{Delta}``                 |
    +                 +                         +---------------------------------------------------+
    |                 |                         | ``Msg{Delta}`` and ``Msg{Gaussian}``              |
    +                 +                         +---------------------------------------------------+
    |                 |                         | ``Msg{Delta}`` and ``Msg{Gamma}``                 |
    +                 +                         +---------------------------------------------------+
    |                 |                         | ``Msg{Delta}`` and ``Msg{InvGamma}``              |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ``Msg{Beta}``           | ``Msg{Beta}`` and ``Msg{Beta}``                   |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ``Msg{Gamma}``          | ``Msg{Gamma}`` and ``Msg{Gamma}``                 |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ``Msg{Gaussian}``       | ``Msg{Gaussian}`` and ``Msg{Gaussian}``           |
    +                 +-------------------------+---------------------------------------------------+
    |                 | ``Msg{Gaussian}`` \*    | ``Msg{Gaussian}`` and ``Msg{StudentsT}``          |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ``Msg{Bernoulli}``      | ``Msg{Bernoulli}`` and ``Msg{Bernoulli}``         |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ``Msg{InvGamma}``       | ``Msg{InvGamma}`` and ``Msg{InvGamma}``           |
    +-----------------+-------------------------+-------------------------+-------------------------+

.. type:: ExponentialNode

    ::

          in        out
        ----->[exp]----->

    :Node function: ``f(in,out) = δ(out - exp(in))``
    :Interfaces:    1. ``i[:in]``, 2. ``i[:out]``
    :Construction:  ``ExponentialNode(id="something")``

    Message computation rules (\* = approximation):

    +-----------------+---------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface            |
    + Rule            +-------------------------+-------------------------+
    |                 | 1                       | 2                       |
    +=================+=========================+=========================+
    | sumProduct!     | ↑↓ ``Msg{Delta}``       | ↑↓ ``Msg{Delta}``       |
    +                 +-------------------------+-------------------------+
    |                 | ↑↓ ``Msg{Gaussian}`` \* | ↑↓ ``Msg{Gamma}`` \*    |
    +-----------------+-------------------------+-------------------------+

.. type:: FixedGainNode

    ::

          in      out
        ----->[A]----->

    :Node function: ``f(in,out) = δ(A*in-out)``
    :Interfaces:    1 ``1[:in]``, 2. ``i[:out]``
    :Construction:  ``FixedGainNode(A::Matrix, id="something")``

    Message computation rules:

    +-----------------+---------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface            |
    + Rule            +-------------------------+-------------------------+
    |                 | 1                       | 2                       |
    +=================+=========================+=========================+
    | sumProduct!     | ↑↓ ``Msg{Delta}``       | ↑↓ ``Msg{Delta}``       |
    +                 +-------------------------+-------------------------+
    |                 | ↑↓ ``Msg{Gaussian}``    | ↑↓ ``Msg{Gaussian}``    |
    +-----------------+-------------------------+-------------------------+


.. type:: GaussianNode

    ::

               mean
                |
                v  out
         ----->[N]----->
        precision/
        variance

    :Node function: ``f(mean,variance,out) = N(out|mean,variance)``
    :Interfaces:    1. ``i[:mean]``, 2. ``i[:variance]`` or ``i[:precision]``, 3. ``i[:out]``
    :Construction:  ``GaussianNode(id="something", form=:moment, m=optional, V=optional)``

    The ``GaussianNode`` outputs a Gaussian distribution from variable mean and variable variance or precision. Upon construction the role of the second interface is set to represent a variance or precision by setting the ``form`` argument to ``:moment or ``:precision`` respectively. The ``m`` and ``V`` arguments allow the user to fix the value for the mean and/or variance interface. Fixed interfaces are not explicitly created.

    Message computation rules:

    +-------------+-------------------------------------------------------------------------+
    |             | Input (↓) and output (↑) per interface                                  |
    + Rule        +---------------------+-----------------------------+---------------------+
    |             | 1                   | 2                           |  3                  |
    +=============+=====================+=============================+=====================+
    | sumProduct! | ↑ ``Msg{Gaussian}`` | ↓ ``Msg{Delta}``            | ↓ ``Msg{Delta}``    |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Msg{Delta}``    | ↑ ``Msg{(Inv)Gamma}``       | ↓ ``Msg{Delta}``    |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Msg{Delta}``    | ↓ ``Msg{Delta}``            | ↑ ``Msg{Gaussian}`` |
    +-------------+---------------------+-----------------------------+---------------------+
    | vmp!        | ↑ ``Msg{Gaussian}`` | ↓ ``(Inv)Gamma``            | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Gaussian``      | ↑ ``Msg{(Inv)Gamma}``       | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Gaussian``      | ↓ ``(Inv)Gamma``            | ↑ ``Msg{Gaussian}`` |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↑ ``Msg{StudentsT}``| ↓ ``Msg{Gamma}``            | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Msg{Gaussian}`` | ↑ ``Msg{Gamma}``            | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``NormalGamma``                                 | ↑ ``Msg{Gaussian}`` |
    +-------------+---------------------+-----------------------------+---------------------+


.. type:: SigmoidNode

    ::

         real     bin
        ----->[σ]----->

    :Node function: ``f(real,bin) = σ(real ⋅ bin)``
    :Interfaces:    1. ``i[:real]``, 2. ``i[:bin]``
    :Construction:  ``SigmoidNode(id="something")``

    The SigmoidNode links a real-valued variable to a binary (∈ {-1,+1}) one.

    Message computation rules:

    +-----------------+---------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface            |
    + Rule            +-------------------------+-------------------------+
    |                 | 1                       | 2                       |
    +=================+=========================+=========================+
    | sumProduct!     |  ↓ ``Msg{Delta}``       | ↑  ``Msg{Bernoulli}``   |
    +                 +-------------------------+-------------------------+
    |                 |  ↓ ``Msg{Gaussian}`` \* | ↑  ``Msg{Bernoulli}``   |
    +-----------------+-------------------------+-------------------------+
    | ep!             | ↑  ``Msg{Gaussian}``    |  ↓ ``Msg{Delta{Bool}}`` |
    +                 +-------------------------+-------------------------+
    |                 | ↑  ``Msg{Gaussian}``    |  ↓ ``Msg{Bernoulli}``   |
    +-----------------+-------------------------+-------------------------+


.. type:: TerminalNode

    (alias ``PriorNode``)
    ::

             out
        [T]----->

    :Node function: ``f(out) = T.value``
    :Interfaces:    1. ``i[:out]``
    :Construction:  ``TerminalNode(value, id="something")``

    A ``TerminalNode`` is used to terminate an edge. It forces the variable represented by the connected edge to ``value``. The terminal node always emits a ``Message`` with payload ``value`` (which is a :class:`ProbabilityDistribution`). It can be used to introduce priors or data into the factor graph.

    Message computation rules:

    +-----------------+---------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface            |
    + Rule            +---------------------------------------------------+
    |                 | 1                                                 |
    +=================+===================================================+
    | sumProduct!     | ↑↓ ``Msg{Any}``                                   |
    +-----------------+-------------------------+-------------------------+


Combined nodes
~~~~~~~~~~~~~~

.. type:: GainAdditionNode

    Combines a :class:`FixedGainNode` with an :class:`AdditionNode` for higher computational efficiency::

                 | in1
                 |
             ____|____
             |   v   |
             |  [A]  |
             |   |   |
         in2 |   v   | out
        -----|->[+]--|---->
             |_______|

    :Node function: ``f(in1,in2,out) = δ(out - A*in1 - in2)``
    :Interfaces:    1. ``i[:in1]``, 2. ``i[:in2]``, 3. ``i[:out]``
    :Construction:  ``GainAdditionNode(A, id="something")``

    Message computation rules:

    +-----------------+-----------------------------------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface                                      |
    + Rule            +-------------------------+-------------------------+-------------------------+
    |                 | 1                       | 2                       |  3                      |
    +=================+=========================+=========================+=========================+
    | sumProduct!     | ↓↑ ``Msg{Gaussian}``    | ↓↑ ``Msg{Gaussian}``    | ↓↑ ``Msg{Gaussian}``    |
    +-----------------+-------------------------+-------------------------+-------------------------+

.. type:: GainEqualityNode

    Combines a :class:`FixedGainNode` with an :class:`EqualityNode` for higher computational efficiency::

             _________
         in1 |       | in2
        -----|->[=]<-|-----
             |   |   |
             |   v   |
             |  [A]  |
             |___|___|
                 | out
                 v

    :Node function: ``f(in1,in2,out) = δ(in1 - A*out)*δ(in2 - A*out)``
    :Interfaces:    1. ``i[:in1]``, 2. ``i[:in2]``, 3. ``i[:out]``
    :Construction:  ``GainEqualityNode(A, id="something")``

    Message computation rules:

    +-----------------+-----------------------------------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface                                      |
    + Rule            +-------------------------+-------------------------+-------------------------+
    |                 | 1                       | 2                       |  3                      |
    +=================+=========================+=========================+=========================+
    | sumProduct!     | ↓↑ ``Msg{Gaussian}``    | ↓↑ ``Msg{Gaussian}``    | ↓↑ ``Msg{Gaussian}``    |
    +-----------------+-------------------------+-------------------------+-------------------------+

.. _composite-nodes:

Composite nodes
---------------

It is possible to create a node that contains an internal :class:`FactorGraph` to define the node function. Such a node is called a ``CompositeNode``.

.. type:: CompositeNode

    A ``CompositeNode`` behaves like a normal ``Node`` from the outside, but contains an *internal graph* that defines the node function. The interfaces of a ``CompositeNode`` are linked to :class:`TerminalNode` instances in its internal graph. A ``CompositeNode`` can easily be constructed from a :class:`FactorGraph`, and it allows one to build hierarchical models since the internal graph may contain other composite nodes.
    ::

        type CompositeNode <: Node
            id::Symbol
            interfaces::Array{Interface,1}
            i::Dict{Symbol,Interface}
            internal_graph::FactorGraph
            # ... and some internal stuff
        end

    Field ``i`` contains named interface handles, for example ``comp_node.i[:out]`` might be identical to ``comp_node.interfaces[2]``. To create multiple instances of a ``CompositeNode``, use ``copy(node::Node, id_of_copy::Symbol)``.



Wrapping a FactorGraph in a CompositeNode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most straightforward way of constructing a ``CompositeNode`` is to first build its internal factor graph, and then wrapping this graph in a ``CompositeNode``. This can be achieved using the following constructor::

    CompositeNode(graph::FactorGraph=current_graph, terminals...; id=generateNodeId(), deterministic=false)

Here, ``terminals`` is an array of :class:`TerminalNode` instances in ``graph`` that should be linked to interfaces of the created ``CompositeNode``. The name of a linked ``TerminalNode`` determines the name of the corresponding :class:`Interface`. Once the ``graph`` is wrapped in a newly created ``CompositeNode``, a new ``FactorGraph`` is created which contains the new ``CompositeNode`` as its only node. This new ``FactorGraph`` becomes the current graph. Example::

    # Build CompositeNode with node function f(in,out) = δ(out - 3*in)

    # Step 1: build internal graph
    g = FactorGraph()
    TerminalNode(id=:t_in)
    TerminalNode(3.0, id=:t_constant)
    TerminalNode(id=:t_out)
    AdditionNode(id=:adder)
    Edge(n(:t_in), n(:adder).i[:in1])
    Edge(n(:t_constant), n(:adder).i[:in2])
    Edge(n(:adder).i[:out], n(:t_out))

    # Step 2: wrap graph in CompositeNode, link t_in & t_out to interfaces
    CompositeNode(g, t_in, t_out, id=:comp_add3) # Creates a new FactorGraph that contains the constructed CompositeNode

    # Step 3: build higher-level graph
    TerminalNode(id=:in)
    TerminalNode(id=:out)
    Edge(n(:t_in), n(:comp_add3).i[:in])
    Edge(n(:comp_add3).i[:out], n(:t_out))

When more instances of the composite node are needed, the ``copy`` function may come in handy::

    copy(n(:comp_add_3), id=:comp_add_3_copy)

will create a new composite node instance with the (optional) id ``:comp_add_3_copy`` that can immediately be used for further construction.


Message computation rules
~~~~~~~~~~~~~~~~~~~~~~~~~

Since a ``CompositeNode`` behaves like a normal ``Node`` from the outside, one can just call a message calculation rule like :func:`sumProduct!` on it. The message will in general be calculated by performing message passing on the internal graph of the composite node. If no suitable custom calculation rule is defined in the ``CompositeNode``, ForneyLab will try to automatically derive a suitable :class:`Algorithm` on the internal graph to calculate the desired message. However, this might not be possible or desireable, for example if the internal graph contains loops. In such cases it is required to define a *custom message calculation rule* using the function ``addRule!()``.

.. function:: addRule!(composite_node::CompositeNode, outbound_interface::Interface, message_calculation_rule::Function, algorithm::Algorithm)

    Add a custom message calculation rule to ``composite_node``. The outbound interface for the rule is specified by ``outbound_interface``. The type of message calculation rule is specified in ``message_calculation_rule``, and can be any valid message calculation rule, like :func:`sumProduct!` or :func:`vmp!`. The ``algorithm`` argument contains the :class:`Algorithm` that yields the desired :class:`Message`.

Note that it's possible to define so called *shortcut rules* using ``addRule!()``. One might for example implement the ``sumProduct!`` rule for a specific interface of a ``CompositeNode`` as a closed-form equation that is derived by hand instead of performing sum-product message passing on the internal factor graph.


.. seealso::
    **Demo:** `Composite nodes <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/06_composite_nodes.ipynb>`_
