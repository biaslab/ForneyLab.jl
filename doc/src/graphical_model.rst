*******************************
 Setting up the graphical model
*******************************

.. seealso::
    **Demo:** `Basics <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/01_basics.ipynb>`_

Building the graphical model amounts to adding nodes and edges to a :class:`FactorGraph`. This chapter introduces the main types and functions that are involved in setting up the graphical model. For some examples, have a look at the following demos: ``random_walk``, ``simple_kalman``.

The ``FactorGraph`` type
========================

.. type:: FactorGraph

    A directed graph consisting of factor nodes and edges::

        type FactorGraph
            n::Dict{Symbol, Node} # Nodes
            e::Dict{Symbol, Edge} # Edges
            # and some internal stuff for attaching buffers etc...
        end

    Upon loading ForneyLab, an empty `FactorGraph` is constructed.
    To create a new one, use ``FactorGraph()``. 
    There is always one *globally active* ``FactorGraph`` instance. A newly constructed :class:`Node` or :class:`Edge` is always added to the current active :class:`FactorGraph`. 

The following functions are available to get/set the currently active factor graph:

.. function:: currentGraph()

    Returns the currently active :class:`FactorGraph` instance.

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

The ``id`` fields holds a unique id, which can be passed to the constructor as a keyword argument. The ``interfaces`` field holds a list of :class:`Interface` objects. An :class:`Edge` connects two interfaces of different nodes to each other. The ``i`` field stores named handles to the interfaces, i.e. ``gain_node.i[:out]`` is equivalent to ``gain_node.interfaces[2]``.

The calling signature of a node constructor depends on the specific type of the node, e.g., 

    AdditionNode(id=:my_adder)  # Node func.: out = in1 + in2
    FixedGainNode(3.0, id=:times_3) # Node func.: out = 3.0 * in1

Nodes in the current graph can be accessed through the function ``node(id::Symbol)`` (which is aliased by the function ``n(id::Symbol)``), e.g.::     

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
            marginal::Union(ProbabilityDistribution, Nothing)
            distribution_type::DataType
        end

    An edge represents a variable, so the ``marginal`` field may contain the marginal :class:`ProbabilityDistribution` over that variable. The ``distribution_type`` field indicates the allowed distribution type of the variable. 

    In general, an ``Edge`` is constructed by passing the tail and head interfaces as well as the distribution type::

        Edge(n(:node1).i[:out], n(:node2).i[:in], GammaDistribution, id=:my_edge)

    If the distribution type is omitted, a :class:`GaussianDistribution` is assumed. For nodes that only have one interface (i.e. :class:`TerminalNode`) or that are symmetrical (i.e. :class:`EqualityNode`), it is also possible to pass the node instead of the interface, e.g., 

        Edge(TerminalNode(), EqualityNode())

    In such cases the constructor will automatically pick the first free interface of the node.

Strictly speaking, a factor graph edge does not need to be directed. However, in ForneyLab all edges are directed to have a consistent meaning for terms like "forward message", "backward messages", and "forward pass". Apart from that, the edge direction has no functional consequences.

ForneyLab does not allow half-edges: every :class:`Edge` should be connected to two nodes at all times. Open ended edges should be terminated by a :class:`TerminalNode`. 

Edges in the current graph can be accessed through the function ``edge(id::Symbol)`` (which is aliased by the function ``e(id::Symbol)``), e.g.::

    edge(:my_edge)
    e(:my_edge)


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

Note that the `Wrap` function only takes *terminal* nodes as arguments.  

.. type:: Wrap(source::TerminalNode, sink::TerminalNode, id::Symbol)

    Constructs a wrap from ``source`` to ``sink`` in the currently active graph and can optionally be given an ``id``.

.. function:: wrap(id::Symbol)

    Returns the ``wrap`` with the corresponding ``id``.

.. function:: wraps(graph::FactorGraph)

    Returns a set of all wraps present in ``graph``. If ``graph`` is omitted, the currently active graph is assumed.

.. function:: wraps(node::TerminalNode)

    Returns a set of all wraps in which ``node`` is involved. Note that a node can be the source in multiple wraps, but it can be a sink at most once.

.. function:: clearWraps(graph::FactorGraph)

    Remove all wraps from :class:`FactorGraph` ``graph``. If ``graph`` is omitted, the currently active graph is assumed.


Interfacing to and from the graph
=================================

.. seealso::
    **Demo:** `Simple Kalman filter <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/04_simple_kalman.ipynb>`_

There are several helper functions that enable the user to connect the graph with the outside world. Reading input and writing output is done through buffers. Several helper functions are available to reset buffers and messages in the graph.

Input to the graph
------------------

Read buffers hold input data that is read into the graph from the outside world. The data is stored in a ``buffer`` vector that is coupled with a terminal ``node``. Upon each call of the :func:`step()` function, the first element of each read buffer is moved to the value field of their coupled nodes.

.. function:: setReadBuffer(node::TerminalNode, buffer::Vector)

    Couples the vector ``buffer`` as read buffer to the :class:`TerminalNode` ``node``.

.. function:: setReadBuffer(nodes::Vector{TerminalNode}, buffer::Vector)

    Couples a read buffer to a batch of nodes. This function can be used to couple input data with a graph that models multiple (time) slices, such as a (mini-)batch. Upon each :func:`step()` call, a number of elements equal to the length of the ``nodes`` vector is moved from the beginning of ``buffer`` to the ``nodes`` value fields (in their respective order). 

Output from the graph
---------------------

Write buffers push message payloads and marginals on a specific interface or edge to an output vector. Upon definition, these functions return an empty output buffer that grows upon each call to :func:`step()`.

.. function:: buffer = setWriteBuffer(interface::Interface)

    Pushes the message payload on ``interface`` to ``buffer`` upon each step.

.. function:: buffer = setWriteBuffer(edge::Edge)

    Pushes the marginal distribution on ``edge`` to ``buffer`` upon each step.

Resetting the graph
-------------------

.. function:: clearBuffers()

    Removes all couplings with read and write buffers.

.. function:: emptyWriteBuffers()

    Resets all write buffers to an empty vector. Pointers to the write buffers are preserved.

.. function:: clearMessages!()

    Clears all messages in the graph.
