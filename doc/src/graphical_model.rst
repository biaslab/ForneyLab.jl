*******************************
 Setting up the graphical model
*******************************

.. seealso::
    **Demo:** `Basics <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/01_basics.ipynb>`_

Building the graphical model amounts to adding nodes and edges to a :class:`FactorGraph`. This chapter introduces the main types and functions that are involved in setting up the graphical model.

The ``FactorGraph`` type
========================

.. type:: FactorGraph

    A directed graph consisting of factor nodes and edges::

        type FactorGraph
            nodes::Dict{Symbol, Node} # Nodes
            edges::Dict{Symbol, Edge} # Edges
            # and some internal stuff for attaching buffers etc...
        end

The following functions are available to get/set the currently active ``FactorGraph``:

.. function:: currentGraph()

    Returns the currently active :class:`FactorGraph` instance. A new ``FactorGraph`` is constructed if there is none.

.. function:: setCurrentGraph(graph::FactorGraph)

    Make ``graph`` the currently active factor graph.


Nodes
=====

A node in a :class:`FactorGraph` is always of a subtype of ``abstract Node``. ForneyLab comes with a bunch of built-in nodes for commonly used functions, such as the :class:`AdditionNode` and the :class:`EqualityNode`. The :doc:`nodes` chapter describes these in more detail, as well as how to build custom nodes. It's important to know that a node should contain a couple of required fields::

    type MinimalNode <: Node
        id::Symbol
        interfaces::Array{Interface,1}
        i::Dict{Symbol, Interface}
    end

The ``id`` field holds a unique id, which can be passed to the constructor as a keyword argument. The ``interfaces`` field holds a list of :class:`Interface` objects. An :class:`Edge` connects two interfaces of different nodes to each other. The ``i`` field stores named handles to the interfaces, i.e. ``gain_node.i[:out]`` is equivalent to ``gain_node.interfaces[2]``.

The calling signature of a node constructor depends on the specific type of the node, e.g.::

    AdditionNode(id=:my_adder)  # Node func.: out = in1 + in2
    GainNode(gain=3.0, id=:times_3) # Node func.: out = 3.0 * in1

A ``Node`` constructor always adds the constructed node to the current graph. To delete a ``Node`` from a :class:`FactorGraph`, use ``delete!(graph::FactorGraph, node::Node)``. Nodes in the current graph can be accessed through the function ``node(id::Symbol)`` (which is aliased by the function ``n(id::Symbol)``), e.g.::

    node(:my_adder)
    n(:my_adder)


The ``Edge`` type
=================

.. type:: Edge

    An ``Edge`` is directed and connects two :class:`Interface` instances of different nodes::

        type Edge <: AbstractEdge
            # [tail]------>[head]
            id::Symbol
            tail::Interface
            head::Interface
            marginal::Union(ProbabilityDistribution, Void)
            distribution_type::DataType
        end

    An edge represents a variable, so the ``marginal`` field may contain the marginal :class:`ProbabilityDistribution` over that variable. The ``distribution_type`` field indicates the allowed distribution type of the variable.

    In general, an ``Edge`` is constructed by passing the tail and head interfaces as well as the distribution type::

        Edge(n(:node1).i[:out], n(:node2).i[:in], GammaDistribution, id=:my_edge)

    If the distribution type is omitted, a :class:`GaussianDistribution` is assumed. For nodes that only have one interface (i.e. :class:`TerminalNode`) or that are symmetrical (i.e. :class:`EqualityNode`), it is also possible to pass the node instead of the interface, e.g.,::

        Edge(TerminalNode(), EqualityNode())

    In such cases the constructor will automatically pick the first free interface of the node.

    The ``Edge`` constructor will add the edge to the current graph (the head and tail nodes should already belong to that graph). To delete an ``Edge`` from a :class:`FactorGraph`, use ``delete!(graph::FactorGraph, edge::Edge)``:

    .. function:: delete!(graph::FactorGraph, edge::Edge)

        Delete the specified ``Edge`` from ``graph``.


Strictly speaking, a factor graph edge does not need to be directed. However, in ForneyLab all edges are directed to have a consistent meaning for terms like "forward message", "backward messages", and "forward pass". Apart from that, the edge direction has no functional consequences.

ForneyLab does not allow half-edges: every :class:`Edge` should be connected to two nodes at all times. Open ended edges should be terminated by a :class:`TerminalNode`.

Edges in the current graph can be accessed through the function ``edge(id::Symbol)`` (which is aliased by the function ``eg(id::Symbol)``), e.g.::

    edge(:my_edge)
    eg(:my_edge)


Example
=======

Consider the following simple factor graph::

          | C1    | C2
          |       |
      X1  v   X2  v   X3
    ---->[+]---->[+]---->

ForneyLab does not allow 'half-edges' that are connected to just one node. Instead, half-edges should be terminated by a :class:`TerminalNode`. Taking this into account, one could implement this factor graph as follows::

    g = FactorGraph()

    # Create nodes
    TerminalNode(id=:t_x1)
    TerminalNode(id=:t_c1)
    TerminalNode(id=:t_c2)
    TerminalNode(id=:t_x3)
    AdditionNode(id=:adder_1)
    AdditionNode(id=:adder_2)

    # Create edges
    Edge(n(:t_x1), n(:adder_1).i[:in1])
    Edge(n(:t_c1), n(:adder_1).i[:in2])
    Edge(n(:adder_1).i[:out], n(:adder_2).i[:in1])
    Edge(n(:t_c2), n(:adder_2).i[:in2])
    Edge(n(:adder_2).i[:out], n(:t_x3))

Chaining factor graph sections
==============================

.. seealso::
    **Demo:** `Random walk <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/03_random_walk.ipynb>`_

In practical situations it is common for a factor graph to be a concatination of identical sections. In such cases it might not be necessary to build the entire factor graph explictly. Instead, it is possible to just build one section, and define how the sections are chained together. This can be done in ForneyLab by defining *wraps*::

    # Random walk chain
    #          | C
    #          |
    #    X[n]  v  X[n+1]
    # ...---->[+]-------> ...

    g = FactorGraph()
    TerminalNode(id=:X_prev)
    TerminalNode(id=:X_next)
    TerminalNode(id=:C)
    AdditionNode(id=:adder)

    Edge(n(:X_prev), n(:adder).[:in1])
    Edge(n(:C), n(:adder).[:in2])
    Edge(n(:adder).i[:out], n(:X_next))

    Wrap(n(:X_next), n(:X_prev)) # X_next becomes X_prev in the next section


.. type:: Wrap

    A ``Wrap`` is a special kind of :class:`Edge` that connects to :class:`TerminalNode` instances such that the involved :class:`FactorGraph` is 'folded'. Inbound messages towards the *source* terminal of a ``Wrap`` will be tranferred to the *sink* of that ``Wrap`` by the :func:`step` function.

    .. function:: Wrap(source::TerminalNode, sink::TerminalNode; id::Symbol)

        Constructs a wrap from ``source`` to ``sink`` in the currently active graph. The keyword argument ``id`` is optional.

    .. function:: wrap(id::Symbol, graph::FactorGraph=currentGraph())

        Returns the ``Wrap`` in ``graph`` with the specified ``id``.

    .. function:: wraps(graph::FactorGraph, graph::FactorGraph=currentGraph())

        Returns a set of all ``Wrap`` instances present in ``graph``.

    .. function:: wraps(node::TerminalNode, graph::FactorGraph=currentGraph())

        Returns a set of all ``Wrap`` instances in which ``node`` is involved. Note that a node can be the source in multiple wraps, but it can be a sink at most once.

    .. function:: delete!(graph::FactorGraph, wrap::Wrap)

        Delete the specified ``Wrap`` from ``graph``.


Interfacing to and from the graph
=================================

.. seealso::
    **Demo:** `Kalman filter <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/04_simple_kalman.ipynb>`_

To link a :class:`FactorGraph` to the outside world, so-called buffers can be used. A buffer can be used to insert data into the graph ('read buffer') or to extract data from the graph ('write buffer'). Some helper functions are available to work with these buffers.

Input to the graph
------------------

Read buffers hold input data that is read into the graph from the outside world. The data is stored in a ``buffer`` vector that is attached to a :class:`TerminalNode`. Every time the :func:`step()` function is called, the first element of each read buffer is moved to the value field of the corresponding terminal node. The following functions are available to attach and detach read buffers:

.. function:: attachReadBuffer(node::TerminalNode, buffer::Vector)

    Attaches the vector ``buffer`` as a read buffer to the :class:`TerminalNode` ``node``.

.. function:: attachReadBuffer(nodes::Vector{TerminalNode}, buffer::Vector)

    Attaches a read buffer to a batch of nodes. This function can be used to couple input data with a graph that models multiple (time) slices, such as a (mini-)batch. On each call to :func:`step()`, a number of elements equal to the length of the ``nodes`` vector is moved from the beginning of ``buffer`` to the value fields of ``nodes`` (in their respective order).

.. function:: detachReadBuffer(node::TerminalNode)

    Detach the read buffer from ``node``.

Output from the graph
---------------------

Write buffers allow message payloads and edge marginals to be extracted from the :class:`FactorGraph`. A write buffer is a ``Vector{ProbabilityDistribution}``, and can be attached to either an :class:`Interface` or an :class:`Edge`. Every call to :func:`step()` will result in exactly one element (message payload or marginal) being pushed onto every write buffer. The following functions are available:

.. function:: attachWriteBuffer(interface::Interface)

    Returns an empty write buffer attached to ``interface``. Every call to :func:`step` will result in the payload of the outbound message on ``interface`` being pushed to the buffer.

.. function:: detachWriteBuffer(interface::Interface)

    Detaches the write buffer attached to ``interface``.

.. function:: attachWriteBuffer(edge::Edge)

    Returns an empty write buffer attached to ``edge``. Every call to :func:`step` will result in the marginal on ``edge`` being pushed to the buffer.

.. function:: detachWriteBuffer(edge::Edge)

    Detaches the write buffer attached to ``edge``.

Resetting the graph
-------------------

.. function:: detachBuffers(graph::FactorGraph=currentGraph())

    Detaches all read and write buffers.

.. function:: emptyWriteBuffers(graph::FactorGraph=currentGraph())

    Truncates the contents of all write buffers.

.. function:: clearMessages!(graph::FactorGraph=currentGraph())

    Clears all messages in the graph.
