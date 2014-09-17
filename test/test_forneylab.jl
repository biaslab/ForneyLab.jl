# This file contains the general ForneyLab tests.
# Tests for specific node and message types are
# found in test_nodes.jl and test_messages.jl 

module TestForneyLab

using FactCheck
using ForneyLab

# Test style and helpers
include("test_style.jl") # Test style conventions on source files
include("integration_helpers.jl") # Helper file for integration tests, contains backgrounds and validations
include("test_helpers.jl") # Tests for ForneyLab helper methods

# Probability distribution test
facts("General ProbabilityDistribution unit tests") do
    for probdist_type in subtypes(ProbabilityDistribution)
        context("$(probdist_type) should have a default constructor and a == operator") do
            @fact probdist_type()==probdist_type() => true
        end
    end
end

# Distribution tests
include("distributions/test_gaussian.jl")
include("distributions/test_gamma.jl")
include("distributions/test_inverse_gamma.jl")
include("distributions/test_normal_gamma.jl")
include("distributions/test_students_t.jl")

# Basic building blocks and methods tests
include("test_node.jl")
include("test_edge.jl")
include("test_graph.jl")
include("test_factorization.jl")

# Nodes
include("nodes/test_addition.jl")
include("nodes/test_terminal.jl")
include("nodes/test_equality.jl")
include("nodes/test_fixed_gain.jl")
include("nodes/test_gaussian.jl")

# Composite nodes
include("nodes/composite/test_gain_addition.jl")
include("nodes/composite/test_gain_equality.jl")
include("nodes/composite/test_linear.jl")

# Marginal calculations
include("test_calculate_marginal.jl")

# Message passing
include("test_message_passing.jl")
include("test_generate_schedule.jl")

# (S)VMP tests
include("test_vmp.jl")

try
    # Try to load user-defined extensions tests
    include("$(Main.FORNEYLAB_EXTENSION_DIR)/test/test_forneylab_extensions.jl")
end

end # module TestForneyLab