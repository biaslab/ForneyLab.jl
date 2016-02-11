module ForneyLab

using Optim
using YAML
using LaTeXStrings

export ProbabilityDistribution, UnivariateProbabilityDistribution, MultivariateProbabilityDistribution
export sumProductRule!, expectationRule!, variationalRule!
export vague, self, ==, isProper, sample, dimensions
export setVerbosity
export prepare!

# Verbosity
verbose = false
setVerbosity(is_verbose=true) = global verbose = is_verbose
printVerbose(msg) = if verbose println(msg) end

# ForneyLab helpers
include("helpers.jl")

# Other includes
import Base.show, Base.convert

# High level abstracts
abstract AbstractEdge # An Interface belongs to an Edge, but Interface is defined before Edge. Because you can not belong to something undefined, Edge will inherit from AbstractEdge, solving this problem.
abstract ProbabilityDistribution # ProbabilityDistribution can be carried by a Message or an Edge (as marginal)
abstract UnivariateProbabilityDistribution <: ProbabilityDistribution
abstract MultivariateProbabilityDistribution <: ProbabilityDistribution

# Node
include("node.jl")

# Message type
include("message.jl")

# Extract dimensionality from message or distribution (note exception for normal-gamma in normal_gamma.jl)
dimensions{T<:MultivariateProbabilityDistribution}(message::Message{T}) = typeof(message.payload).parameters[end]
dimensions(distribution::MultivariateProbabilityDistribution) = typeof(distribution).parameters[end]
dimensions{T<:MultivariateProbabilityDistribution}(message_type::Type{Message{T}}) = message_type.parameters[1].parameters[end]
dimensions{T<:MultivariateProbabilityDistribution}(distribution_type::Type{T}) = distribution_type.parameters[end]

dimensions{T<:UnivariateProbabilityDistribution}(message::Message{T}) = 1
dimensions(distribution::UnivariateProbabilityDistribution) = 1
dimensions{T<:UnivariateProbabilityDistribution}(message_type::Type{Message{T}}) = 1
dimensions{T<:UnivariateProbabilityDistribution}(distribution_type::Type{T}) = 1


# Univariate distributions
include("distributions/univariate/delta.jl")
include("distributions/univariate/bernoulli.jl")
include("distributions/univariate/gaussian.jl")
include("distributions/univariate/gamma.jl")
include("distributions/univariate/inverse_gamma.jl")
include("distributions/univariate/students_t.jl")
include("distributions/univariate/beta.jl")
include("distributions/univariate/log_normal.jl")

# Multivariate distributions
include("distributions/multivariate/mv_delta.jl")
include("distributions/multivariate/mv_gaussian.jl")
include("distributions/multivariate/normal_gamma.jl")
include("distributions/multivariate/wishart.jl")
include("distributions/multivariate/mv_log_normal.jl")

# Basic ForneyLab building blocks and methods
include("interface.jl")
include("edge.jl")
include("schedule.jl")

# Nodes
include("nodes/addition.jl")
include("nodes/terminal.jl")
include("nodes/equality.jl")
include("nodes/gain.jl")
include("nodes/gaussian.jl")
include("nodes/exponential.jl")
include("nodes/gain_addition.jl")
include("nodes/gain_equality.jl")
include("nodes/sigmoid.jl")

# Graph, wraps and algorithm
include("factor_graph.jl")
include("wrap.jl")
include("inference_algorithm.jl")

# Composite nodes
include("nodes/composite.jl")

# Methods for calculating marginals
include("distributions/calculate_marginal.jl")

# Generic methods
include("message_passing.jl")
include("step.jl")

# Utils
include("visualization.jl")

# InferenceAlgorithms
include("algorithms/sum_product/sum_product.jl")
include("algorithms/loopy_sum_product/loopy_sum_product.jl")
include("algorithms/variational_bayes/variational_bayes.jl")
include("algorithms/expectation_propagation/expectation_propagation.jl")

# Shared preparation methods for inference algorithms
include("algorithms/preparation.jl")

vague{T<:UnivariateProbabilityDistribution}(dist_type::Type{T}) = vague!(T())

include("docstrings.jl")

end # module ForneyLab
