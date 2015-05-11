**************
 Algorithms
**************

.. type:: Algorithm


Execution
=========

.. function:: execute(algorithm, graph)


.. function:: step(algorithm, graph)


.. function:: run(algorithm, graph)


Schedules
=========

.. type:: Schedule


.. type:: ScheduleEntry


The sum-product algorithm
=========================

ForneyLab implements inference algorithms in separate modules. Most importantly, the `SumProduct` module holds specific `SumProduct.Algorithm` constructors and an automatic scheduler for generating a sumproduct message passing schedule.

Upon algorithm construction the generated schedule is stored under the `:schedule` key of the `algorithm.fields` dictionary.

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


Variational message passing
===========================

The variational message passing (VMP) module implements VMP as described by Dauwels in his 2007 paper "On variational message passing on factor graphs". The module is capable of conducting both mean field and structured VMP and implements several algorithm specific constructors, an auto scheduler and several helper types required for execution.

The q-factorization is stored under the `:factorization` key of the `algorithm.fields` dictionary and references the different subgraphs. The actual q-distributions are stored under the `:q_distributions` key and the number of iterations under `:n_iterations`.

Algorithm constructors
----------------------

.. function:: VMP.Algorithm(::FactorGraph)

    Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on the argument graph, with a mean field factorization.

.. function:: VMP.Algorithm(cluster_edges...)

    Generates a vmp algorithm to calculate the messages towards write buffers and timewraps defined on the argument graph, with a structured factorization. When unpacked, (the extension of) the elements of type `Set{Edge}` in `cluster_edges` define the separate subgraphs.

An optional field `n_iterations=50` specifies the number of VMP iterations.

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

    The `QFactorization` type stores the variational factorization of the graph. The `edge_to_subgraph` attribute contains a dictionary for fast subgraph lookup::

        type QFactorization
            factors::Array{Subgraph, 1}
            edge_to_subgraph::Dict{Edge, Subgraph}
        end

.. type:: QDistribution

    The `QDistribution` type stores local q-distributions that are the approximate marginals on the external edges. The `edges` attribute defined the set of edges on which `distribution` is defined::

        type QDistribution
            distribution::ProbabilityDistribution
            edges::Set{Edge} # Edges on which the distribution is defined
        end

