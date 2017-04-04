module ForneyLab

# using Optim, LaTeXStrings

# export ProbabilityDistribution, Univariate, Multivariate, MatrixVariate, AbstractDelta
# export sumProductRule!, expectationRule!, variationalRule!
# export InferenceAlgorithm
# export vague, self, ==, isProper, sample, dimensions, pdf, logpdf, mean, var, cov
export setVerbosity
# export prepare!
# export rules

# Verbosity
verbose = false
setVerbosity(is_verbose=true) = global verbose = is_verbose
printVerbose(msg) = if verbose println(msg) end

# Helpers
include("helpers.jl")
include("dependency_graph.jl")

# Other includes
import Base.show, Base.convert, Base.==, Base.mean, Base.var, Base.cov

# High level abstracts
abstract AbstractEdge # An Interface belongs to an Edge, so AbstractEdge has to be defined before Interface
abstract AbstractVariable
abstract InferenceAlgorithm

# Low-level internals
include("factor_node.jl")
include("interface.jl")
include("probability_distribution.jl")
include("approximation.jl")
include("edge.jl")
include("variable.jl")

# Factor nodes
include("factor_nodes/constant.jl")
include("factor_nodes/equality.jl")
include("factor_nodes/gaussian.jl")
include("factor_nodes/gaussian_mean_variance.jl")

# include("nodes/equality.jl")
# include("nodes/addition.jl")
# include("nodes/terminal.jl")
# include("nodes/gain.jl")
# include("nodes/gaussian.jl")
# include("nodes/gaussian_mixture.jl")
# include("nodes/exponential.jl")
# include("nodes/gain_addition.jl")
# include("nodes/gain_equality.jl")
# include("nodes/sigmoid.jl")
# include("nodes/bernoulli.jl")
# include("nodes/categorical.jl")



# Graph and algorithm
include("factor_graph.jl")
# include("inference_algorithm.jl")

# Composite nodes
include("factor_nodes/composite.jl")

# Generic methods
include("message_passing.jl")
# include("step.jl")

# Update rules
include("update_rules/constant.jl")
include("update_rules/gaussian_mean_variance.jl")

# Utils
include("visualization.jl")

# InferenceAlgorithms
include("algorithms/sum_product/sum_product.jl")
# include("algorithms/loopy_sum_product/loopy_sum_product.jl")
# include("algorithms/variational_bayes/variational_bayes.jl")
# include("algorithms/expectation_propagation/expectation_propagation.jl")

# # Shared preparation methods for inference algorithms
# include("algorithms/preparation.jl")

# vague{T<:Univariate}(dist_type::Type{T}) = vague!(T())
# *(x::ProbabilityDistribution, y::ProbabilityDistribution) = prod!(x, y) # * operator for probability distributions
# *(x::Message, y::Message) = prod!(x.payload, y.payload) # * operator for messages

# include("docstrings.jl")

end # module ForneyLab
