using ForneyLab

# Initial settings
N1              = 50                                           # Number of observed samples first cluster
N2              = 50                                           # Number of observed samples second cluster
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

#Nodes
GaussianMixtureNodePar(id=:gm)
EqualityNode(id=:m1_eq)
EqualityNode(id=:w1_eq)
EqualityNode(id=:pi_eq)
TerminalNode(Delta(),id=:y)
PriorNode(Bernoulli(0.5),id=:z)
PriorNode(PartitionedDistribution([Gaussian(m=50.0,V=12.0),Gaussian(m=6.0,V=12.0)]),id=:m1_min_1)
PriorNode(PartitionedDistribution([Gamma(a=1,b=0.01),Gamma(a=1, b=0.01)]),id=:w1_min_1)
PriorNode(Beta(a=2.0,b=2.0),id=:pi_min_1)
TerminalNode(vague(PartitionedDistribution{Gaussian, 2}),id=:m1_n)
TerminalNode(vague(PartitionedDistribution{Gamma,2}),id=:w1_n)
TerminalNode(vague(Beta),id=:pi_n)

#Edges
Edge(n(:pi_eq).i[1],n(:gm).i[:pi],id=:pi)
Edge(n(:m1_eq).i[1],n(:gm).i[:m],id=:m1)
Edge(n(:w1_eq).i[1],n(:gm).i[:w],id=:w1)
Edge(n(:z).i[:out],n(:gm).i[:z],id=:z)
Edge(n(:y).i[:out],n(:gm).i[:x],id=:y)
Edge(n(:m1_min_1).i[:out], n(:m1_eq).i[:2])
Edge(n(:w1_min_1).i[:out], n(:w1_eq).i[:2])
Edge(n(:pi_min_1).i[:out], n(:pi_eq).i[:2])
Edge(n(:m1_n).i[:out], n(:m1_eq).i[:3])
Edge(n(:w1_n).i[:out], n(:w1_eq).i[:3])
Edge(n(:pi_n).i[:out], n(:pi_eq).i[:3])

Wrap(n(:pi_n),n(:pi_min_1))
Wrap(n(:m1_n),n(:m1_min_1))
Wrap(n(:w1_n),n(:w1_min_1))

#attach the observed data
attachReadBuffer(n(:y), deepcopy(y));

# attach write buffers to the wanted variables
m1_est = attachWriteBuffer(n(:m1_eq).i[3])
w1_est = attachWriteBuffer(n(:w1_eq).i[3])
pi_est = attachWriteBuffer(n(:pi_eq).i[3]);

# Specify the variational algorithm for n_its vmp iterations
algo = VariationalBayes(Dict(   eg(:m1) => PartitionedDistribution{Gaussian,2},
                                eg(:w1) => PartitionedDistribution{Gamma,2},
                                eg(:z)  => Bernoulli,
                                eg(:pi) => Beta,
                                eg(:y)  => Gaussian),
                        n_iterations=n_its)

show(algo)

#run the algorithm
run(algo)

#print the true variables
# println("True mean 1: $(true_mean1)")
# println("True precision 1: $(1/true_variance1)")
# println("True mean 2: $(true_mean2)")
# println("True precision 2: $(1/true_variance2)")
# println("Number of samples: $(N1+N2)")
# println("Sample mean 1: $(round(mean(y_observations1),2))")
# println("Sample precision 1: $(round(1/var(y_observations1),2))")
# println("Sample mean 2: $(round(mean(y_observations2),2))")
# println("Sample precision 2: $(round(1/var(y_observations2),2))")
# println("Total sample mean: $(round(mean(y),2))")
# println("Total sample precision: $(round(1/var(y),2))")
#
# #print the estimated variables
# println("\n----- Online estimation after $(n_its) VMP updates per sample -----")
println("Mean estimate 1: $(round(mean(m1_est[end].factors[1])[1],2)), with variance $(round(var(m1_est[end].factors[1])[1,1],2))")
println("Precision estimate 1: $(round(mean(w1_est[end].factors[1]),2))")
println("Mean estimate 2: $(round(mean(m1_est[end].factors[2])[1],2)), with variance $(round(var(m1_est[end].factors[2])[1,1],2))")
println("Precision estimate 2: $(round(mean(w1_est[end].factors[2]),2))")
