**************
 Algorithms
**************

Algorithms that are to be executed on a :class:`FactorGraph`, are instances of the ``InferenceAlgorithm`` type. ForneyLab implements ``SumProduct``, ``LoopySumProduct``, ``VariationalBayes`` and ``ExpectationPropagation`` algorithms. While algorithm specific details differ, the general mechanisms for construction and execution remain the same. For example, the ``SumProduct`` algorithm is defined as follows.

.. type:: SumProduct

    A specific algorithm type for sum-product message passing::

        type SumProduct
            execute::Function
            schedule::Schedule
        end

    The ``SumProduct`` algorithm holds a ``schedule``, which is executed upon calling the ``execute`` function.

ForneyLab implements several constructors for algorithm types. Construction of an algorithm generally consists of these steps:

1. Generate a schedule towards specific interfaces, wraps and/or write buffers;
2. Infer the message distribution types for each entry in the schedule.

Inference of the distribution types in step 2 allows for efficient execution of the update rules later on.


Execution
=========

After construction, an algorithm can be executed. For efficiency reasons, computations are in-place operations performed on the messages stored on the interfaces. Instead of making a new ``ProbabilityDistribution`` instance, the parameters of the existing ``message.payload`` are altered. Algorithm execution generally consists of the following steps:

1. Ensure that initial outbound messages are present on each outbound interface in the schedule;
2. For each entry in the schedule, pre-compile the update rule;
3. In turn execute each entry in schedule, in-place updating the pre-set messages.

Steps 1 and 2 are the preparation phase. The pre-compilation in step 2 ensures that functions ``run`` and ``step`` only need to call a single pre-compiled update function for each entry in the schedule.

The following functions are exported:

.. function:: execute(algorithm)

    Executes ``algorithm`` on the current graph.


.. function:: step(algorithm, direction::Symbol)

    Buffer and wrap aware version of ``execute(algorithm)``. Performs the following:

    1. Read an element from all read-buffers defined in the current graph
    2. ``execute(algorithm)``
    3. Write an element to all write-buffers defined in the current graph
    4. Propagate messages along all wraps defined in current graph in the corresponding direction

    Before calling ``step``, the algorithm must first be prepared by calling ``prepare!(algorithm)``. Preparing the algorithm ensures that initial messages are set and all update rules are pre-compiled.
    Shorthand version of the function ``step(algorithm)`` exists and is an alias for ``step(algorithm, :forward)``

.. function:: run(algorithm; direction::Symbol)

    Prepares the algorithm and calls ``step(algorithm)`` repeatedly until one of the read-buffers is exhausted or until the end of the block is reached. Keyword argument ``direction`` is optional and by default takes ``:forward`` value.


Schedules
=========

Algorithms executed on factor graphs are usually based on *message passing*. In message passing, a sequence of messages is calculated, where the next message usually depends on previously calculated messages. A ``Schedule`` defines a sequence of message calculations.

.. type:: Schedule

    Defines a sequence of :class:`Message` calculations::

        typealias Schedule Vector{ScheduleEntry}

    .. function:: execute(schedule::Schedule)

        Executes all entries in ``schedule`` once.


.. type:: ScheduleEntry

    Specifies a message calculation operation::

        type ScheduleEntry
            node::Node
            outbound_interface_id::Int64
            rule::Function  # Refers to the general message calculation rule; for example sumProductRule! or variationalRule!.
            execute::Function # Compiled rule call: () -> rule(node, Val{outbound_interface_id}, rule_arguments...).

            # And some omitted fields
        end

    The ``ScheduleEntry`` is the workhorse of ForneyLab. Most importantly, the ``execute`` field holds the pre-compiled (anonymous) function for the message update. All other fields are simply there to facilitate the proper construction of ``execute``. The ``execute`` function is called upon execution of the ``ScheduleEntry``.


The sum-product algorithm
=========================

.. seealso::
    **Demo:** `Kalman filter <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/04_simple_kalman.ipynb>`_

The ``SumProduct`` algorithm comes with several constructors and an automatic scheduler for generating a sum-product message passing schedule.

Algorithm constructors for sum-product message passing only work for acyclic graphs. For graphs with cycles, the ``LoopySumProduct`` algorithm can be used.

.. function:: SumProduct(::FactorGraph)

    Generates a sum-product algorithm with a schedule towards interfaces connected to write buffers and wraps.

.. function:: SumProduct(::Interface)

    Defines a sum-product algorithm with a schedule towards the argument interface.

.. function:: SumProduct(::Vector{Interface})

    Generates a sum-product algorithm that at least propagates messages to all interfaces in the argument vector.

.. function:: SumProduct(::Edge)

    Defines a sum-product algorithm with a schedule towards the forward and backward interfaces of the argument edge and calculates the corresponding marginal.


Automatic scheduler
-------------------

.. function:: generateSumProductSchedule(::FactorGraph)

    Returns a sum-product message passing schedule that passes messages towards interfaces concerning write buffers and wraps as defined by the argument graph. The scheduler works through depth first search.


The loopy sum-product algorithm
===============================

.. seealso::
    **Demo:** `Loopy belief propagation <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/05_loopy_belief_propagation.ipynb>`_

The ``LoopySumProduct`` algorithm is similar to the ``SumProduct`` algorithm, but then for graphs with cycles.

.. function:: LoopySumProduct(::FactorGraph; breaker_messages=Dict{Interface, Message}(), n_iterations=50, ...)

    Constructs a loopy sum-product algorithm that propagates to defined write buffers and wraps. Breaker messages specified by the ``breaker_messages`` dictionary are pre-set on the corresponding interfaces. From there a standard sum-product message passing schedule is generated. Upon execution, this schedule is repeated for ``n_iterations``.

.. function:: LoopySumProduct(::Interface; breaker_messages=Dict{Interface, Message}(), n_iterations=50, ...)

    Constructs a loopy sum-product algorithm towards an interface.


The variational message passing algorithm
=========================================

.. seealso::
    | **Demo:** `Naive variational message passing <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/07_naive_variational_message_passing.ipynb>`_
    | **Demo:** `Structured variational message passing <https://github.com/spsbrats/ForneyLab.jl/blob/master/demo/08_structured_variational_message_passing.ipynb>`_

The ``VariationalBayes`` algorithm implements variational message passing (VMP) as described by Dauwels in his 2007 paper "On variational message passing on factor graphs". The algorithm supports both mean field and structured VMP. ForneyLab implements several algorithm specific constructors, an auto scheduler and several helper types required for execution.

The factorization of the recognition distribution is stored under the ``factorization`` field of the algorithm and references the different subgraphs. The actual recognition distributions are stored under the ``recognition_distributions`` field and the number of iterations under ``:n_iterations``.


Algorithm constructors
----------------------

.. function:: VariationalBayes(recognition_distribution_types::Dict, ::FactorGraph; n_iterations=50)

    Generates a VMP algorithm to calculate the messages towards write buffers and timewraps defined on the argument graph, with a as specified by the ``recognition_distribution_types`` dictionary.

The factorization of the recognition distribution is specified by the edge(array)-to-distribution-type dictionary called ``recognition_distribution_types``. The conventions for passing the recognition distribution factorization are best specified by example.

The snippet below specifies a full (mean field) factorization around a Gaussian node::

    algo = VariationalBayes(Dict(
        eg(:mean) => GaussianDistribution,
        eg(:prec) => GammaDistribution,
        eg(:out)  => GaussianDistribution))

For defining a full factorization over multiple graph sections, edges with similar distributions are grouped in a column vector::

    algo = VariationalBayes(Dict(
        [eg(:mean1), eg(:mean2)] => GaussianDistribution,
        [eg(:prec1), eg(:prec2)] => GammaDistribution,
        [eg(:out1),  eg(:out2) ] => GaussianDistribution))

Edges belonging to the same cluster are grouped in the rows of a matrix. The following snippet specifies a joint recognition distribution over the mean and precision (note the lack of a separating comma)::

    algo = VariationalBayes(Dict(
        [eg(:mean) eg(:prec)] => NormalGammaDistribution,
         eg(:out)             => GaussianDistribution))

For more examples, consult the VMP demos.


Automatic scheduler
-------------------

.. function:: generateVariationalBayesSchedule!(::RecognitionFactorization, ::FactorGraph)

    Generates and stores an (internal and external) schedule for VMP.


VMP specific types
------------------

.. type:: Subgraph

    The internal edges of subgraphs are non-overlapping clusters, which together define the factorization of the recognition distribution. The VMP algorithm executes updates for the subgraphs (corresponding with the factors) in turn::

        type Subgraph
            internal_edges::Set{Edge}
            internal_schedule::Schedule # Schedule for internal message passing
            external_schedule::Array{Node, 1} # Schedule for marginal updates
        end

.. type:: RecognitionFactorization

    The ``RecognitionFactorization`` type stores the variational factorization of the graph. The ``edge_to_subgraph`` attribute contains a dictionary for fast subgraph lookup::

        type RecognitionFactorization
            factors::Array{Subgraph, 1}
            edge_to_subgraph::Dict{Edge, Subgraph}
        end

.. type:: RecognitionDistribution

    The ``RecognitionDistribution`` type stores local recognition distributions. The ``edges`` attribute defined the local set of edges on which ``distribution`` is defined::

        type RecognitionDistribution
            distribution::ProbabilityDistribution
            edges::Set{Edge} # Edges on which the distribution is defined
        end


The expectation propagation algorithm
=====================================

The ``ExpectationPropagation`` algorithm automatically derives an expectation propagation message passing algorithm. The expectation propagation (EP) algorithm is similar to (loopy) belief propagation as implemented by the sum-product algorithm. For some nodes, the exact sum-product messages cannot be expressed analytically in the desired form, rendering the sum-product algorithm unusable. In these cases, the EP algorithm provides a solution by projecting the 'difficult' messages on the family of desired distributions. The interfaces that generate the 'difficult' messages are called sites. The outbound messages on the sites are called "expectations", and represent local approximations to the 'true' messages. The inbound messages on the sites are called "cavity distributions", and they capture the effect of the rest of the graph (usually prior + other sites) on the marginal. Since the expectation message depends on the cavity distribution, the EP algorithm creates implicit loops in the factor graph. Because of this, the EP message passing schedule has to be executed multiple times for the messages to converge.

The expectation messages on the sites are calculated by the :func:`expectationRule!`. This rule should be implemented for all nodes connected to sites. In contrast to :func:`sumProductRule!`, :func:`expectationRule!` also consumes the inbound message on the outbound interface (site).


.. function:: ExpectationPropagation(sites::Vector{Interface}; ...)

    Generates an EP algorithm to incrementally approximate the marginal distributions of the variables (edges) connected to the specified 'sites'. The generated message passing schedule will respect the order of the sites.
    The following optional keyword arguments may be passed:

    - ``num_iterations``: a positive integer indicating the maximum number of iterations (default=100).
    - ``callback``: a function that is called after each iteration. This function can be used for example to check converge or to collect intermediate results. If the callback function returns ``true``, the algorithm is terminated.
