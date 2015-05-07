*******************************
 Setting up the graphical model
*******************************

Building the graphical model amounts to adding nodes and edges to a :class:`FactorGraph`. This chapter introduces the main types and functions that are involved in setting up the graphical model. For some examples, have a look at the following demos: ``random_walk``, ``simple_kalman``.

The ``FactorGraph`` type
------------------------

.. type:: FactorGraph

    A directed graph consisting of factor nodes and edges::

        type FactorGraph
            nodes::Set{Node}
            edges::Set{Edge}
            # and some internal stuff for attaching buffers etc...
        end

    Upon loading ForneyLab, an empty `FactorGraph` is constructed.
    To create a new one, use ``FactorGraph()``. 
    There is always one *globally active* ``FactorGraph`` instance. An :class:`Edge` is always added to the currently active :class:`FactorGraph` upon construction, along with the nodes connected to it. 

The following functions are available to get/set the currently active factor graph:

.. function:: currentGraph()

    Returns the currently active :class:`FactorGraph` instance.

.. function:: setCurrentGraph(graph::FactorGraph)

    Make ``graph`` the currently active factor graph.


Nodes
-----

A node in a :class:`FactorGraph` is always of a subtype of ``abstract Node``. ForneyLab comes with a bunch of built-in nodes for commonly used functions, such as the :class:`AdditionNode` and the :class:`EqualityNode`. The :doc:`nodes` chapter describes these in more detail, as well as how to build custom nodes. It's important to know that a node should contain a couple of required fields::

    type MinimalNode <: Node
        name::ASCIIString
        interfaces::Array{Interface,1}
    end

The ``name`` fields holds a unique name, which can be passed to the constructor as a keyword argument. The ``interfaces`` field holds a list of :class:`Interface` objects. An :class:`Edge` connects two interfaces of different nodes to eachother. The calling signature of a node constructor depends on the specific type of the node::

    addition_node = AdditionNode(name="my_adder")  # Node func.: out = in1 + in2
    gain_node = FixedGainNode(3.0, name="times_3") # Node func.: out = 3.0 * in1

Some node types provide named interface handles for convenience. I.e. ``gain_node.out`` is equivalent to ``gain_node.interfaces[2]``.

The ``Edge`` type
-----------------

.. type:: Edge

    An ``Edge`` is directed and connects two :class:`Interface` instances of different nodes::

        type Edge <: AbstractEdge
            # [tail]------>[head]
            tail::Interface
            head::Interface
            marginal::Union(ProbabilityDistribution, Nothing)
            distribution_type::DataType
        end

    An edge represents a variable, so the ``marginal`` field may contain the marginal :class:`ProbabilityDistribution` over that variable. The ``distribution_type`` field indicates the allowed distribution type of the variable. 

    In general, an ``Edge`` is constructed by passing the tail and head interfaces as well as the distribution type::

        edge = Edge(node1.out, node2.interfaces[1], GammaDistribution)

    If the distribution type is omitted, a :class:`GaussianDistribution` is assumed. For nodes that only have one interface (i.e. :class:`TerminalNode`) or that are symmetrical (i.e. :class:`EqualityNode`), it is also possible to pass the node instead of the interface::

        edge = Edge(TerminalNode(), EqualityNode())

    In such cases the constructor will automatically pick the first free interface of the node.

Strictly speaking, a factor graph edge does not need to be directed. However, in ForneyLab all edges are directed to have a consistent meaning for terms like "forward message", "backward messages", and "forward pass". Apart from that, the edge direction has no functional consequences.

ForneyLab does not allow half-edges: an :class:`Edge` should be connected to two nodes at all times. Open ended edges should be terminated by a :class:`TerminalNode`. 

Example
-------

Consider the following simple factor graph::

          | C1    | C2           
          |       |       
      X1  v   X2  v   X3 
    ---->[+]---->[+]---->

ForneyLab does not allow 'half-edges' that are connected to just one node. Instead, half-edges should be terminated by a :class:`TerminalNode`. Taking this into account, one could implement this factor graph like::

    g = FactorGraph()

    # Create nodes
    t_x1 = TerminalNode()
    t_c1 = TerminalNode()
    t_c2 = TerminalNode()
    t_x3 = TerminalNode()
    adder_1 = AdditionNode(name="adder_1")
    adder_2 = AdditionNode(name="adder_2")

    # Create edges
    Edge(t_x1, adder_1.in1)
    Edge(t_c1, adder_1.in2)
    Edge(adder_1.out, adder_2.in1)
    Edge(t_c2, adder_2.in2)
    Edge(adder_2.out, t_x3)

Chaining factor graph sections
------------------------------

In practical situations it is common for a factor graph to be a concatination of identical sections. In such cases it might not be necessary to build the entire factor graph explictly. Instead, it is possible to just build one section, and define how the sections are chained together. This can be done in ForneyLab by defining *wraps*::

    # Random walk chain
    #          | C          
    #          |           
    #    X[n]  v  X[n+1]
    # ...---->[+]-------> ...

    g = FactorGraph()
    X_prev = TerminalNode()
    X_next = TerminalNode()
    C = TerminalNode()
    adder = AdditionNode()

    Edge(X_prev, adder.in1)
    Edge(C, adder.in2)
    Edge(adder.out, X_next)

    wrap(X_next, X_prev) # X_next becomes X_prev in the next section


.. function:: wrap(from, to, graph)

    Creates a wrap from :class:`TerminalNode` ``from`` to :class:`TerminalNode` ``to`` in :class:`FactorGraph` ``graph``. If ``graph`` is omitted, the currently active graph is assumed.

.. function:: clearWraps(graph)

    Remove all wraps from :class:`FactorGraph` ``graph``. If ``graph`` is omitted, the currently active graph is assumed.

