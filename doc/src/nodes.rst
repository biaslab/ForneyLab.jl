**************
 Nodes
**************

This chapter describes the implementation of nodes in ForneyLab. For more details, have a look at the source files, such as ``src/nodes/addition.jl``.

.. _node-anatomy:

The anatomy of nodes
--------------------

Each node type is a subtype of abstract ``Node``. A node should contain at least the fields ``name``, ``interfaces`` and ``i``. Let's look at the definition of the built-in :class:`AdditionNode`::

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

A factor node captures a specific *node function*, which involves the variables that are represented by the connected edges. The :class:`AdditionNode` for example captures the addition function: ``f(in1,in2,out) = δ(in1+in2-out)``. When running a message passing algorithm on a :class:`FactorGraph`, the node function specifies how the outbound messages are calculated from the inbound messages. An outbound message is calculated according to a message calculation rule, also simply called *rule*. Rules are implemented for specific nodes and message types, making heavy use of Julia's multiple dispatch system.

ForneyLab supports the following rules:

.. function:: sumProductRule!(node::Node, outbound_interface_index::Type{Val{i}}, outbound_dist::ProbabilityDistribution, inbound_messages...)

    Calculates the outbound message from the incoming messages on the other interfaces according to the sum-product algorithm.
    Example implementation::

        function sumProductRule!(node::AdditionNode,
                                 outbound_interface_index::Type{Val{3}},
                                 outbound_dist::GaussianDistribution,
                                 msg_in1::Message{GaussianDistribution},
                                 msg_in2::Message{GaussianDistribution},
                                 msg_out::Any)
            # Calculate outbound message on interface 3 (:out)
            # according to the sum-product rule

            # Perform computations as in in-place operation on outbound_dist

            return outbound_dist
        end

    The calling signature consists of:

    1. The node;
    2. The index (index in node.interfaces) of the outbound interface as a value type;
    3. The distribution currently present on the outbound interface;
    4. The inbound messages on *all* interfaces of the node (ordered by interface id).

.. function:: variationalRule!(node::Node, outbound_interface_index::Type{Val{i}}, outbound_dist::ProbabilityDistribution, marginals_and_messages...)

    Similar to :func:`sumProductRule!`, but on some interfaces marginals are required instead of messages. This rule is used for variational message passing (vmp).
    Example implementation::

        function variationalRule!(node::GaussianNode,
                                  outbound_interface_index::Type{Val{1}},
                                  outbound_dist::GaussianDistribution,
                                  marg_mean::Any,
                                  marg_prec::GammaDistribution,
                                  marg_y::GaussianDistribution)
            # Calculate outbound message on interface 1 (:mean)
            # according to the variational rule

            # Perform computations as in in-place operation on outbound_dist

            return outbound_dist
        end

    The calling signature consists of:

    1. The node;
    2. The index (index in node.interfaces) of the outbound interface as a value type;
    3. The distribution currently present on the outbound interface;
    4. The inbound messages/marginals on *all* interfaces of the node (ordered by interface id).

.. function:: expectationRule!(node::Node, outbound_interface_index::Type{Val{i}}, outbound_dist::GaussianDistribution, inbound_messages...)

    Similar to :func:`sumProductRule!`, but also the inbound message on the outbound interface is consumed (this messages carries the cavity distrubution). This calculation rule is used in the expectation propagation algorithm.

    The calling signature consists of:

    1. The node;
    2. The index (index in node.interfaces) of the outbound interface as a value type;
    3. The distribution currently present on the outbound interface;
    4. The inbound messages on *all* interfaces of the node (ordered by interface id).

Not all message calculation rules have to be implemented for every node, just the ones that will be used. Similarly, the message calculation rule does not have to be implemented for a specific outbound interface of a specific node if that outbound message never has to be calculated.

.. _built-in-nodes:

Built-in nodes
--------------

The following built-in 'elementary' nodes are available in ForneyLab: :class:`AdditionNode`, :class:`EqualityNode`, :class:`ExponentialNode`, :class:`GainNode`, :class:`GaussianNode`, :class:`SigmoidNode`, :class:`TerminalNode`.

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

    Message computation rules:

    +-----------------+---------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface            |
    + Rule            +-------------------------+-------------------------+
    |                 | 1                       | 2                       |
    +=================+=========================+=========================+
    | sumProduct!     | ↑↓ ``Msg{Delta}``       | ↑↓ ``Msg{Delta}``       |
    +                 +-------------------------+-------------------------+
    |                 | ↑↓ ``Msg{Gaussian}``    | ↑↓ ``Msg{LogNormal}``   |
    +-----------------+-------------------------+-------------------------+

.. type:: GainNode

    ::

            gain
              |
         in   V   out
       ----->[A]----->


    :Node function: ``f(in,out,gain) = δ(out - gain*in)``, where ``gain`` is either provided upon construction of the node and is a fixed value or is supplied via gain interface.
    :Interfaces:    1 ``i[:in]``, 2. ``i[:out]``, 3. ``i[:gain]``
    :Construction:  ``GainNode(gain=[2.0], id="something")`` or ``GainNode(id="something")``

    The ``GainNode`` implements the the multiplication of a random variable with a non-random variable (encoded by a ``DeltaDistribution``).

    Message computation rules:

    +-----------------+-----------------------------------------------------------------------------+
    |                 | Input (↓) and output (↑) per interface                                      |
    + Rule            +-------------------------+-------------------------+-------------------------+
    |                 | 1                       | 2                       | 3.                      |
    +=================+=========================+=========================+=========================+
    | sumProduct!     | ↑↓ ``Msg{Delta}``       | ↑↓ ``Msg{Delta}``       |                         |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↑↓ ``Msg{Gaussian}``    | ↑↓ ``Msg{Gaussian}``    |                         |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↑  ``Msg{Gaussian}``    | ↑  ``Msg{Gaussian}``    | ↑  ``Msg{Delta}``       |
    +                 +-------------------------+-------------------------+-------------------------+
    |                 | ↓  ``Msg{Gaussian}``    | ↓  ``Msg{Gaussian}``    | ↑  ``Msg{Delta}``       |
    +-----------------+-------------------------+-------------------------+-------------------------+


.. type:: GaussianNode

    ::

               mean
                |
                v  out
         ----->[N]----->
        precision/
        log-precision/
        variance/

    :Node function: ``f(mean,variance,out) = N(out|mean,variance)``
    :Interfaces:    1. ``i[:mean]``, 2. ``i[:variance]``, ``i[:log_variance]``, or ``i[:precision]``, 3. ``i[:out]``
    :Construction:  ``GaussianNode(id="something", form=:moment, m=optional, V=optional)``

    The ``GaussianNode`` outputs a Gaussian distribution from variable mean and variable variance or precision. Upon construction the role of the second interface is set to represent a variance, precision or log-precision by setting the ``form`` argument to ``:moment``, ``:precision`` or ``:log_variance`` respectively. The ``m`` and ``V`` arguments allow the user to fix the value for the mean and/or variance interface. Fixed interfaces are not explicitly created.

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
    |             | ↑ ``Msg{Gaussian}`` | ↓ ``Gaussian``              | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Gaussian``      | ↑ ``Msg{(Inv)Gamma}``       | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Gaussian``      | ↑ ``Msg{Gaussian}``         | ↓ ``Gaussian``      |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Gaussian``      | ↓ ``(Inv)Gamma``            | ↑ ``Msg{Gaussian}`` |
    +             +---------------------+-----------------------------+---------------------+
    |             | ↓ ``Gaussian``      | ↓ ``Gaussian``              | ↑ ``Msg{Gaussian}`` |
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

    Combines a :class:`GainNode` with an :class:`AdditionNode` for higher computational efficiency::

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

    Combines a :class:`GainNode` with an :class:`EqualityNode` for higher computational efficiency::

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

A node in which multiple node functions are combined into one node function is called a ``CompositeNode`` In ForneyLab 0.4 the composite node functionality is removed. It will be reintroduced in ForneyLab 0.5.
