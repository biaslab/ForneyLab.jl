module ForneyLab

using Optim
using YAML
using LaTeXStrings

export ProbabilityDistribution, UnivariateProbabilityDistribution, MultivariateProbabilityDistribution
export sumProduct!, ep!, vmp!
export vague, self, ==, isProper, sample
export setVerbosity

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
include("algorithms/variational_bayes/variational_bayes.jl")
include("algorithms/expectation_propagation/expectation_propagation.jl")

# Functions for message post-processing
vague(dist::ProbabilityDistribution) = vague(typeof(dist))

function __init__()
    # Module-global variable to keep track of currently active InferenceAlgorithm
    global current_algorithm = nothing
end

end # module ForneyLab
