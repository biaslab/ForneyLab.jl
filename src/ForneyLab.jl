module ForneyLab

using Optim
using YAML
using LaTeXStrings

export Node, CompositeNode, ProbabilityDistribution
export sumProduct!
export vague, expectation, self, ==
export current_graph

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
show(io::IO, nodes::Union(Set{Node}, Vector{Node})) = [show(io, node) for node in nodes]
abstract CompositeNode <: Node

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

# Composite nodes
include("nodes/composite/gain_addition.jl")
include("nodes/composite/gain_equality.jl")
include("nodes/composite/sigmoid.jl")
include("nodes/composite/general.jl")

# Graphs and factorizations
include("graph.jl")

# Methods for calculating marginals
include("distributions/calculate_marginal.jl")

# Generic methods
include("message_passing.jl")
include("generate_schedule.jl")
include("step.jl")

# Utils
include("visualization.jl")

# Functions for message post-processing
vague(dist::ProbabilityDistribution) = vague(typeof(dist))
expectation(dist::ProbabilityDistribution) = DeltaDistribution(mean(dist))
self(x::Any) = x

end # module ForneyLab
