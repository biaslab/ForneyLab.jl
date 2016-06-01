using ForneyLab

# Initial settings
N1              = 200                                           # Number of observed samples first cluster
N2              = 200                                           # Number of observed samples second cluster
n_its           = 50                                           # Number of vmp iterations
true_mean1      = 30.0                                          # Mean first cluster
true_variance1  = 0.1                                           # Variance first cluster
true_mean2      = 1.0                                           # Mean second cluster
true_variance2  = 0.1                                           # Variance second cluster
y_observations1 = sqrt(true_variance1)*randn(N1) + true_mean1   #observarions first cluster
y_observations2 = sqrt(true_variance2)*randn(N2)  + true_mean2  #observations second cluster

permutations = shuffle(collect(1:N1+N2))
y = [y_observations1; y_observations2][permutations]            # mix datapoint of two clusters

# Build graph
for k=1:(N1+N2)
    GaussianMixtureNode(id=:gm*k) # s() for symbol concatenation
    EqualityNode(id=:m1_eq*k)
    EqualityNode(id=:w1_eq*k)
    EqualityNode(id=:pi_eq*k)
    TerminalNode(Delta(y[k]), id=:y*k) # Observed y values are stored in terminal node values
    PriorNode(ForneyLab.Bernoulli(0.5),id=:z*k)
    Edge(n(:pi_eq*k).i[1],n(:gm*k).i[:pi],id=:pi_e*k)
    Edge(n(:m1_eq*k).i[1],n(:gm*k).i[:m],id=:m1_e*k)
    Edge(n(:w1_eq*k).i[1],n(:gm*k).i[:w],id=:w1_e*k)
    Edge(n(:z*k).i[:out],n(:gm*k).i[:z],id=:z_e*k)
    Edge(n(:y*k).i[:out],n(:gm*k).i[:x],id=:y_e*k)

    if k > 1 # Connect sections
        Edge(n(:m1_eq*(k-1)).i[2], n(:m1_eq*k).i[3])
        Edge(n(:pi_eq*(k-1)).i[2], n(:pi_eq*k).i[3])
        Edge(n(:w1_eq*(k-1)).i[2], n(:w1_eq*k).i[3])
    end
end


PriorNode(PartitionedDistribution([Gaussian(m=50.0,V=12.0),Gaussian(m=6.0,V=5.0)]),id=:m1_start)
PriorNode(PartitionedDistribution([Gamma(a=1,b=0.01),Gamma(a=1, b=0.01)]),id=:w1_start)
PriorNode(ForneyLab.Beta(a=2.,b=2.),id=:pi_start)

Edge(n(:m1_eq*1).i[3],n(:m1_start).i[:out])
Edge(n(:w1_eq*1).i[3],n(:w1_start).i[:out])
Edge(n(:pi_eq*1).i[3],n(:pi_start).i[:out])

TerminalNode(vague(PartitionedDistribution{Gaussian,2}),id=:m1_end)
TerminalNode(vague(PartitionedDistribution{ForneyLab.Gamma,2}),id=:w1_end)
TerminalNode(vague(ForneyLab.Beta),id=:pi_end)

Edge(n(:m1_eq*(N1+N2)).i[2], n(:m1_end))
Edge(n(:w1_eq*(N1+N2)).i[2], n(:w1_end))
Edge(n(:pi_eq*(N1+N2)).i[2], n(:pi_end))

# attach write buffers to the wanted variables
m1_est = attachWriteBuffer(n(:m1_end).i[:out].partner)
w1_est = attachWriteBuffer(n(:w1_end).i[:out].partner)
pi_est = attachWriteBuffer(n(:pi_end).i[:out].partner);

# Specify the variational algorithm for n_its vmp iterations
msg_types = Dict{Interface,DataType}([n(:gm*i).i[:x] => Gaussian for i=1:(N1+N2)])
algo = VariationalBayes(Dict(   eg(:m1_e*(1:(N1+N2))) => PartitionedDistribution{Gaussian,2},
                                eg(:w1_e*(1:(N1+N2))) => PartitionedDistribution{Gamma,2},
                                eg(:z_e*(1:(N1+N2)))  => Bernoulli,
                                eg(:pi_e*(1:(N1+N2))) => Beta,
                                eg(:y_e*(1:(N1+N2)))  => Gaussian),
                        n_iterations=n_its,
                        message_types=msg_types)

show(algo)

#run the algorithm
run(algo);

# #print the true variables
# println("True mean 1: $(true_mean1)")
# println("True precision 1: $(1/true_variance1)")
# println("True mean 2: $(true_mean2)")
# println("True precision 2: $(1/true_variance2)")
# println("Number of samples: $(N1+N2)")
# println("Sample mean 1: $(round(mean(y_observations1),2)\)")
# println("Sample precision 1: $(round(1/var(y_observations1),2))")
# println("Sample mean 2: $(round(mean(y_observations2),2))")
# println("Sample precision 2: $(round(1/var(y_observations2),2))")
# println("Total sample mean: $(round(mean(y),2))")
# println("Total sample precision: $(round(1/var(y),2))")

#print the estimated variables
println("\n----- Online estimation after $(n_its) VMP updates per sample -----")
println("Mean estimate 1: $(round(mean(m1_est[end].factors[1])[1],2)), with variance $(round(var(m1_est[end].factors[1])[1,1],2))")
println("Precision estimate 1: $(round(mean(w1_est[end].factors[1]),2))")
println("Mean estimate 2: $(round(mean(m1_est[end].factors[2])[1],2)), with variance $(round(var(m1_est[end].factors[2])[1,1],2))")
println("Precision estimate 2: $(round(mean(w1_est[end].factors[2]),2))")
