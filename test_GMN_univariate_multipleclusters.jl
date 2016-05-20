using ForneyLab

# Initial settings
N               = [200; 200; 200]                               # Number of observed samples first cluster
n_its           = 200                                     # Number of vmp iterations
true_mean1      = 30.0                                          # Mean first cluster
true_variance1  = 0.1                                           # Variance first cluster
true_mean2      = 1.0                                           # Mean second cluster
true_variance2  = 0.1                                           # Variance second cluster
true_mean3      = 100.
true_variance3  = 0.1
y_observations1 = sqrt(true_variance1)*randn(N[1]) + true_mean1   #observarions first cluster
y_observations2 = sqrt(true_variance2)*randn(N[2])  + true_mean2  #observations second cluster
y_observations3 = sqrt(true_variance3)*randn(N[3])  + true_mean3  #observations second cluster

permutations = shuffle(collect(1:(sum(N))))
y = [y_observations1; y_observations2; y_observations3][permutations]            # mix datapoint of two clusters

#Overwrite the Gamma prior with variables that are equal or bigger than 1
function ForneyLab.vague!(dist::ForneyLab.Gamma)
    dist.a = 1.
    dist.b = 1.
   return dist
end

ForneyLab.vague(::Type{ForneyLab.Gamma}) = ForneyLab.Gamma(a=1., b=1.)
# Build graph

#Nodes
GaussianMixtureNodePar(id=:gm)
EqualityNode(id=:m1_eq)
EqualityNode(id=:w1_eq)
EqualityNode(id=:pi_eq)
TerminalNode(Delta(),id=:y)
PriorNode(Categorical{length(N)}(1/length(N)*ones(length(N))),id=:z)
PriorNode(PartitionedDistribution([Gaussian(m=25.0,V=12.0),Gaussian(m=10.0,V=12.0),Gaussian(m=125.,V=12.)]),id=:m1_min_1)
PriorNode(PartitionedDistribution([Gamma(a=1,b=1.),Gamma(a=1, b=1.), Gamma(a=1,b=1.)]),id=:w1_min_1)
PriorNode(Dirichlet(20.*ones(length(N))),id=:pi_min_1)
TerminalNode(vague(PartitionedDistribution{Gaussian, length(N)}),id=:m1_n)
TerminalNode(vague(PartitionedDistribution{Gamma,length(N)}),id=:w1_n)
TerminalNode(vague(Dirichlet{length(N)}),id=:pi_n)

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
pi_est = attachWriteBuffer(n(:pi_eq).i[3])
z_est  = attachWriteBuffer(n(:gm).i[:z])

# Specify the variational algorithm for n_its vmp iterations
algo = VariationalBayes(Dict(   eg(:m1) => PartitionedDistribution{Gaussian,length(N)},
                                eg(:w1) => PartitionedDistribution{Gamma,length(N)},
                                eg(:z)  => Categorical{length(N)},
                                eg(:pi) => Dirichlet{length(N)},
                                eg(:y)  => Gaussian),
                        n_iterations=n_its)

show(algo)

#run the algorithm
run(algo);


#print the estimated variables
println("\n----- Online estimation after $(n_its) VMP updates per sample -----")
println("Mean estimate 1: $(round(mean(m1_est[end].factors[1])[1],2)), with variance $(round(var(m1_est[end].factors[1])[1,1],2))")
println("Precision estimate 1: $(round(mean(w1_est[end].factors[1]),2))")
println("Mean estimate 2: $(round(mean(m1_est[end].factors[2])[1],2)), with variance $(round(var(m1_est[end].factors[2])[1,1],2))")
println("Precision estimate 2: $(round(mean(w1_est[end].factors[2]),2))")
println("Mean estimate 3: $(round(mean(m1_est[end].factors[3])[1],2)), with variance $(round(var(m1_est[end].factors[3])[1,1],2))")
println("Precision estimate 3: $(round(mean(w1_est[end].factors[3]),2))")
