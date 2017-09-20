module ForneyLab

# Helpers
include("helpers.jl")
include("dependency_graph.jl")

# Other includes
import Base: show, convert, ==, mean, var, cov, *

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
include("factor_nodes/addition.jl")
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

# Utils
include("visualization.jl")

# InferenceAlgorithms
include("algorithms/sum_product/sum_product.jl")
# include("algorithms/loopy_sum_product/loopy_sum_product.jl")
# include("algorithms/variational_bayes/variational_bayes.jl")
# include("algorithms/expectation_propagation/expectation_propagation.jl")

# # Shared preparation methods for inference algorithms
# include("algorithms/preparation.jl")

# Update rules
include("update_rules/constant.jl")
include("update_rules/equality.jl")
include("update_rules/addition.jl")
include("update_rules/gaussian_mean_variance.jl")

*(x::ProbabilityDistribution, y::ProbabilityDistribution) = prod!(x, y) # * operator for probability distributions

# include("docstrings.jl")

# Engines
include("engines/julia/julia.jl")

end # module ForneyLab
