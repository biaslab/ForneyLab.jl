# Run the tests for ForneyLab
# This file is automatically executed if you call Pkg.test("ForneyLab")

# This file contains the general ForneyLab tests.
# Tests for specific node and message types are
# found in test_nodes.jl and test_messages.jl

module TestForneyLab

using FactCheck
using ForneyLab

import Base.==

# Test style and helpers
include("test_style.jl") # Test style conventions on source files
include("integration_helpers.jl") # Helper file for integration tests, contains backgrounds and validations
include("test_helpers.jl") # Tests for ForneyLab helper methods

# Distribution tests
include("test_distributions.jl")
include("distributions/univariate/test_bernoulli.jl")
include("distributions/univariate/test_delta.jl")
include("distributions/univariate/test_gaussian.jl")
include("distributions/univariate/test_gamma.jl")
include("distributions/univariate/test_inverse_gamma.jl")
include("distributions/univariate/test_students_t.jl")
include("distributions/univariate/test_beta.jl")
include("distributions/univariate/test_log_normal.jl")

include("distributions/multivariate/test_mv_delta.jl")
include("distributions/multivariate/test_mv_gaussian.jl")
include("distributions/multivariate/test_normal_gamma.jl")

# Basic building blocks and methods tests
include("test_interface.jl")
include("test_node.jl")
include("test_edge.jl")
include("test_schedule.jl")
include("test_wrap.jl")

# Top level concepts
include("test_factor_graph.jl")

# Node types
include("nodes/test_addition.jl")
include("nodes/test_terminal.jl")
include("nodes/test_equality.jl")
include("nodes/test_fixed_gain.jl")
include("nodes/test_gaussian.jl")
include("nodes/test_exponential.jl")
include("nodes/test_gain_addition.jl")
include("nodes/test_gain_equality.jl")
include("nodes/test_sigmoid.jl")

# Composite nodes
include("nodes/test_composite.jl")

# Message passing
include("test_message_passing.jl")
include("test_step.jl")

# Algorithms
include("test_algorithm.jl")
include("algorithms/sum_product/test_sum_product.jl")
include("algorithms/vmp/test_vmp.jl")
include("algorithms/expectation_propagation/test_expectation_propagation.jl")

exitstatus()

end # module TestForneyLab
