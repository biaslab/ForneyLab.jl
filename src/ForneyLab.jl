module ForneyLab

using Optim
using YAML
using LaTeXStrings

export Node, ProbabilityDistribution
export sumProduct!, vmp!
export vague, self, ==
export current_graph

# Export algorithm modules
export SumProduct
export VMP

# Verbosity
verbose = false
setVerbosity(is_verbose=true) = global verbose = is_verbose
printVerbose(msg) = if verbose println(msg) end

# ForneyLab helpers
include("helpers.jl")

# Other includes
import Base.show, Base.convert

# Top-level abstracts
abstract AbstractEdge # An Interface belongs to an Edge, but Interface is defined before Edge. Because you can not belong to something undefined, Edge will inherit from AbstractEdge, solving this problem.
abstract ProbabilityDistribution # ProbabilityDistribution can be carried by a Message or an Edge (as marginal)
abstract Node
show(io::IO, node::Node) = println(io, "$(typeof(node)) with name $(node.name)")

# Message type
include("message.jl")

# Distributions
include("distributions/delta.jl")
include("distributions/gaussian.jl")
include("distributions/gamma.jl")
include("distributions/inverse_gamma.jl")
include("distributions/normal_gamma.jl")
include("distributions/students_t.jl")
include("distributions/beta.jl")

# Basic ForneyLab building blocks and methods
include("interface.jl")
include("edge.jl")
include("schedule.jl")

# Nodes
include("nodes/addition.jl")
include("nodes/terminal.jl")
include("nodes/equality.jl")
include("nodes/fixed_gain.jl")
include("nodes/gaussian.jl")
include("nodes/exponential.jl")
include("nodes/gain_addition.jl")
include("nodes/gain_equality.jl")
include("nodes/sigmoid.jl")

# Graph and algorithm
include("factor_graph.jl")
include("algorithm.jl")

# Composite nodes
include("nodes/composite.jl")

# Methods for calculating marginals
include("distributions/calculate_marginal.jl")

# Generic methods
include("message_passing.jl")
include("step.jl")

# Utils
include("visualization.jl")

# Algorithms
include("algorithms/sum_product/sum_product.jl")
include("algorithms/vmp/vmp.jl")

# Functions for message post-processing
vague(dist::ProbabilityDistribution) = vague(typeof(dist))
self(x::Any) = x

end # module ForneyLab
