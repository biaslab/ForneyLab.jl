module ForneyLab

using Optim
using LaTeXStrings

export ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, AbstractDelta
export sumProductRule!, expectationRule!, variationalRule!
export InferenceAlgorithm
export vague, self, ==, isProper, sample, dimensions, pdf
export setVerbosity
export prepare!
export rules

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
# Documentation in docstrings.jl
abstract ProbabilityDistribution # ProbabilityDistribution can be carried by a Message or an Edge (as marginal)
abstract Univariate <: ProbabilityDistribution
abstract Multivariate <: ProbabilityDistribution
abstract MatrixVariate <: ProbabilityDistribution

abstract InferenceAlgorithm

# Low-level internals
include("approximation.jl")     # Types related to approximations
include("node.jl")              # Node type
include("message.jl")           # Message type

# Dimensionality of distributions
dimensions(::Univariate) = 1
dimensions(distribution::Multivariate) = typeof(distribution).parameters[end]

# Dimensionality of distribution types
dimensions{T<:Univariate}(::Type{T}) = 1
dimensions{T<:Multivariate}(distribution_type::Type{T}) = distribution_type.parameters[end]

# Dimensionality of messages and message types
dimensions(message::Message) = dimensions(message.payload)
dimensions(message_type::Type{Message}) = dimensions(message_type.parameters[1])

# Delta distributions
include("distributions/univariate/delta.jl")
include("distributions/multivariate/mv_delta.jl")
include("distributions/matrix_variate/matrix_delta.jl")
typealias AbstractDelta Union{Delta,MvDelta,MatrixDelta}

# Univariate distributions
include("distributions/univariate/bernoulli.jl")
include("distributions/univariate/categorical.jl")
include("distributions/univariate/gaussian.jl")
include("distributions/univariate/gamma.jl")
include("distributions/univariate/inverse_gamma.jl")
include("distributions/univariate/students_t.jl")
include("distributions/univariate/beta.jl")
include("distributions/univariate/log_normal.jl")

# Multivariate distributions
include("distributions/multivariate/dirichlet.jl")
include("distributions/multivariate/mv_gaussian.jl")
include("distributions/multivariate/mv_log_normal.jl")
include("distributions/multivariate/normal_gamma.jl")
include("distributions/multivariate/partitioned.jl")

# Matrix variate distributions
include("distributions/matrix_variate/wishart.jl")

# Special distributions
include("distributions/mixture.jl")

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
include("nodes/gaussian_mixture.jl")
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

vague{T<:Univariate}(dist_type::Type{T}) = vague!(T())
*(x::ProbabilityDistribution, y::ProbabilityDistribution) = prod!(x, y) # * operator for probability distributions
*(x::Message, y::Message) = prod!(x.payload, y.payload) # * operator for messages

include("docstrings.jl")

end # module ForneyLab
