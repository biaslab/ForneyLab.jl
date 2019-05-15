module ForneyLab

# Libraries that are not extended
using PDMats

# Helpers
include("helpers.jl")
include("dependency_graph.jl")

# Other includes
import Base: show, convert, ==, *
import Base.Meta: parse
import Base64: base64encode
import LinearAlgebra: diag, det, tr, cholesky, pinv, Adjoint, Transpose
import SparseArrays: spzeros
import SpecialFunctions: digamma, lgamma, lbeta, erfc
import Statistics: mean, var, cov


# High level abstracts
abstract type AbstractEdge end # An Interface belongs to an Edge, so AbstractEdge has to be defined before Interface
abstract type AbstractVariable end
abstract type InferenceAlgorithm end

# Low-level internals
include("factor_node.jl")
include("interface.jl")
include("probability_distribution.jl")
include("edge.jl")
include("variable.jl")

# Factor nodes
include("factor_nodes/clamp.jl")
include("factor_nodes/equality.jl")
include("factor_nodes/addition.jl")
include("factor_nodes/multiplication.jl")
include("factor_nodes/exponential.jl")
include("factor_nodes/gaussian_mean_variance.jl")
include("factor_nodes/gaussian_mean_precision.jl")
include("factor_nodes/gaussian_weighted_mean_precision.jl")
include("factor_nodes/gaussian.jl")
include("factor_nodes/gamma.jl")
include("factor_nodes/log_normal.jl")
include("factor_nodes/wishart.jl")
include("factor_nodes/bernoulli.jl")
include("factor_nodes/categorical.jl")
include("factor_nodes/contingency.jl")
include("factor_nodes/transition.jl")
include("factor_nodes/beta.jl")
include("factor_nodes/dirichlet.jl")
include("factor_nodes/gaussian_mixture.jl")
include("factor_nodes/sigmoid.jl")
include("factor_nodes/nonlinear.jl")
include("factor_nodes/dot_product.jl")

# Factor graph
include("factor_graph.jl")

# Composite nodes
include("factor_nodes/composite.jl")

# Generic methods
include("message_passing.jl")
include("marginals.jl")

# Utils
include("visualization.jl")

# InferenceAlgorithms
include("algorithms/sum_product/sum_product.jl")
include("algorithms/variational_bayes/recognition_factorization.jl")
include("algorithms/variational_bayes/joint_marginals.jl")
include("algorithms/variational_bayes/naive_variational_bayes.jl")
include("algorithms/variational_bayes/structured_variational_bayes.jl")
include("algorithms/expectation_propagation/expectation_propagation.jl")

# Update rules
include("update_rules/clamp.jl")
include("update_rules/equality.jl")
include("update_rules/addition.jl")
include("update_rules/multiplication.jl")
include("update_rules/exponential.jl")
include("update_rules/gaussian_mean_variance.jl")
include("update_rules/gaussian_mean_precision.jl")
include("update_rules/gamma.jl")
include("update_rules/log_normal.jl")
include("update_rules/wishart.jl")
include("update_rules/bernoulli.jl")
include("update_rules/categorical.jl")
include("update_rules/transition.jl")
include("update_rules/beta.jl")
include("update_rules/dirichlet.jl")
include("update_rules/gaussian_mixture.jl")
include("update_rules/sigmoid.jl")
include("update_rules/nonlinear.jl")
include("update_rules/dot_product.jl")

*(x::ProbabilityDistribution, y::ProbabilityDistribution) = prod!(x, y) # * operator for probability distributions

# include("docstrings.jl")

# Engines
include("engines/julia/julia.jl")

end # module ForneyLab
