**************
 Algorithms
**************

Algorithms that are to be executed on a :class:`FactorGraph`, are instances of the ``Algorithm`` type.

.. type:: Algorithm

    Generic type to hold algorithm-specific code and data structures::

        type Algorithm
            execute::Function
            fields::Dict{Symbol, Any}
        end

    The ``execute`` field holds a function that contains the algorithm's code. Algorithms might require specific data structures, such as message passing schedules. Those are saved in the ``fields`` dictionary.

One can construct an ``Algorithm`` object manually, or use one of the available algorithm constructors to build an ``Algorithm`` automatically. For tree-like factor graphs, it is for example possible automatically derive a sum-product message passing algorithm using::

    algo = SumProduct.Algorithm(graph)

See `built-in algorithm generators`_ for more info.


Execution
=========

An :class:`Algorithm` instance can be executed on a :class:`FactorGraph`. The following functions are available:


.. function:: execute(algorithm, graph)

    Executes ``algorithm`` on ``graph``. 


.. function:: step(algorithm, graph)

    Buffer and wrap aware version of ``execute(algorithm)``. Performs the following:

    1. Read an element from all read-buffers defined in ``graph``
    2. ``execute(algorithm, graph)``
    3. Write an element to all write-buffers defined in ``graph``
    4. Apply all wraps defined in ``graph``

.. function:: run(algorithm, graph)

    Calls ``step(algorithm, graph)`` repeatedly until a read-buffer is exhausted.

    
Schedules
=========

.. seealso::
    **Demo:** `Schedules in loopy graphs <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/06_loopy_graphs.ipynb>`_

Algorithms executed on factor graphs are usually based on *message passing*. In message passing, a sequence of messages is calculated, where the next message usually depends on previously calculated messages. A ``Schedule`` defines a sequence of message calculations.

.. type:: Schedule

    Defines as linear sequence of :class:`Message` calculations::

        typealias Schedule Array{ScheduleEntry, 1}

    .. function:: execute(schedule::Schedule)

        Executes all entries in ``schedule`` once.


.. type:: ScheduleEntry

    Specifies a message calculation operation::

        type ScheduleEntry
            interface::interface                # Calculate outbound message on this interface
            message_calculation_rule::Function  # Called to calculate the message. Default is sumProduct!.
            post_processing::Function           # Optional: a function that performs post-processing on the message.
        end

    When a ``ScheduleEntry`` is executed, the ``message_calculation_rule`` function is called. The ``interface`` field specifies the interface on which the outbound message should be calculated. If ``post_processing`` is defined, the payload of the calculated message is past through it. Examples of commonly used post-processing functions are :func:`sample` and :func:`mean`.

    .. function:: execute(schedule_entry::ScheduleEntry)

        Performs the message calculation specified in ``schedule_entry``. The calculated message is saved on the interface specified by ``schedule_entry`` and is also returned. If a post-processing function is defined, it is applied to the result. If the output of the post-processing function is not a :class:`ProbabilityDistribution`, a :class"`DeltaDistribution` is constructed to hold the output.


Built-in algorithm generators
=============================

ForneyLab includes commonly used inference algorithms, which are implemented in submodules. Currently, the following algorithm submodules are available: :ref:`sumproduct-submodule`, :ref:`vmp-submodule`.

.. _sumproduct-submodule:

The sum-product algorithm
=========================

.. seealso::
    **Demo:** `Simple Kalman filter <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/04_simple_kalman.ipynb>`_

The ``SumProduct`` submodule implements the sum-product algorithm. Most importantly, the submodule holds specific ``SumProduct.Algorithm`` constructors and an automatic scheduler for generating a sum-product message passing schedule.

Upon algorithm construction the generated schedule is stored in the ``:schedule`` key of the ``algorithm.fields`` dictionary.

Algorithm constructors
----------------------

Algorithm constructors for sum product message passing only work for acyclic graphs, or for cyclic graphs with pre-set breaker messages.

.. function:: SumProduct.Algorithm(::FactorGraph)

    Generates a sumproduct algorithm with a schedule towards interfaces concerning write buffers and wraps as defined by the argument graph.

.. function:: SumProduct.Algorithm(::Interface)

    Defines a sumproduct algorithm with a schedule towards the argument interface.

.. function:: SumProduct.Algorithm(::Vector{Interface})

    Generates a sumproduct algorithm that at least propagates messages to all interfaces in the argument vector.

.. function:: SumProduct.Algorithm(::Edge)

    Defines a sumproduct algorithm with a schedule towards the forward and backward interfaces of the argument edge and calculates the corresponding marginal.

Automatic scheduler
-------------------

.. function:: SumProduct.generateSchedule(::FactorGraph)

    Returns a sum product message passing schedule that passes messages towards interfaces concerning write buffers and wraps as defined by the argument graph. The scheduler works through depth first search and terminates when it encounters an interface that carries a message. Normally the automatic scheduler can only works for acyclic graphs, so before schedule generation cycles should be broken with breaker messages.  


.. _vmp-submodule:

Variational message passing
===========================

.. seealso::
    | **Demo:** `Naive variational message passing <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/07_naive_variational_message_passing.ipynb>`_
    | **Demo:** `Structured variational message passing <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/08_structured_variational_message_passing.ipynb>`_

The variational message passing (VMP) submodule implements VMP as described by Dauwels in his 2007 paper "On variational message passing on factor graphs". The module is capable of conducting both mean field and structured VMP and implements several algorithm specific constructors, an auto scheduler and several helper types required for execution.

The q-factorization is stored under the ``:factorization`` key of the ``algorithm.fields`` dictionary and references the different subgraphs. The actual q-distributions are stored under the ``:q_distributions`` key and the number of iterations under ``:n_iterations``.

Algorithm constructors
----------------------

.. function:: VMP.Algorithm(::FactorGraph)

    Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on the argument graph, with a mean field factorization.

.. function:: VMP.Algorithm(cluster_edges...)

    Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on the argument graph, with a structured factorization. When unpacked, (the extension of) the elements of type ``Set{Edge}`` in ``cluster_edges`` define the separate subgraphs.

An optional field ``n_iterations=50`` specifies the number of VMP iterations.

Automatic scheduler
-------------------

.. function:: VMP.generateSchedule!(::Subgraph)

    Generates and stores an (internal and external) schedule for VMP on the argument subgraph. Messages within a subgraph are propagated towards wraps, write buffers and external edges. 

VMP specific types
------------------

.. type:: Subgraph

    The internal edges of subgraphs are non-overlapping clusters, which together define the q-factorization. The VMP algorithm executes updates for the subgraphs in turn::

        type Subgraph
            internal_edges::Set{Edge}
            internal_schedule::Schedule # Schedule for internal message passing
            external_schedule::Array{Node, 1} # Schedule for marginal updates
        end

.. type:: QFactorization

    The ``QFactorization`` type stores the variational factorization of the graph. The ``edge_to_subgraph`` attribute contains a dictionary for fast subgraph lookup::

        type QFactorization
            factors::Array{Subgraph, 1}
            edge_to_subgraph::Dict{Edge, Subgraph}
        end

.. type:: QDistribution

    The ``QDistribution`` type stores local q-distributions that are the approximate marginals on the external edges. The ``edges`` attribute defined the set of edges on which ``distribution`` is defined::

        type QDistribution
            distribution::ProbabilityDistribution
            edges::Set{Edge} # Edges on which the distribution is defined
        end

