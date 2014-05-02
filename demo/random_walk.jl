############################################
# Random walk demo
############################################
# Description:
#   Models a random walk with Gaussian noise.
#
#      [N]         [N]
#       |           |
#       v  X(n-1)   v  X(n)
#   -->[+]-------->[+]------>
#
############################################

module RandomWalkDemo

using ForneyLab

n_steps = 5 # Length of the series

noise_nodes = ConstantNode[] # Initialize an empty set of noise nodes [N]
addition_nodes = AdditionNode[] # Initialize an empty set of addition nodes [+]
for step = 1:n_steps
    # Initialize nodes
    push!(noise_nodes, ConstantNode(GaussianMessage(m=[0.0], V=[1.0]), name="noise_node_$(step)")) # Create a Gaussian noise node
    push!(addition_nodes, AdditionNode(name="addition_node_$(step)")) # Create an addition node

    # Connect edges
    Edge(noise_nodes[step].out, addition_nodes[step].in2) # Connect the noise on the addition input
    if step > 1
        Edge(addition_nodes[step-1].out, addition_nodes[step].in1) # Connect the previous addition node to the current
    end
end

# Initial message
initial_node = ConstantNode(GaussianMessage(m=[0.0], V=[0.0])) # The initial message is observed
Edge(initial_node.out, addition_nodes[1].in1) # Connect the initial node to the beginning of the series

msg_out = calculateMessage!(addition_nodes[n_steps].out) # Calculate the message on the final edge
println("The outgoing message on the final edge is a $(typeof(msg_out)), with mean $(msg_out.m) and variance $(msg_out.V)")

end # module