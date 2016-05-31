using ForneyLab
# Initial settings
N1              = 200 # Number of observed samples per cluster
N2             = 200
n_its          = 100 # Number of vmp iterations
true_mean1      = 30.0
true_variance1  = 0.1
true_mean2     = 1.0
true_variance2 = 0.1
y_observations1 = sqrt(true_variance1)*randn(N1) + true_mean1 # y observation buffer
y_observations2 = sqrt(true_variance2)*randn(N2)  + true_mean2

permutations = shuffle(collect(1:N1+N2))
y = [y_observations1; y_observations2][permutations]

# Build graph
GaussianMixtureNodePar(id=:gm)
EqualityNode(id=:m1_eq)
#EqualityNode(id=:m2_eq)
EqualityNode(id=:w1_eq)
#EqualityNode(id=:w2_eq)
EqualityNode(id=:pi_eq)
TerminalNode(Delta(),id=:y)
PriorNode(Bernoulli(0.5),id=:z)
PriorNode(PartitionedDistribution([Gaussian(m=27.0,V=12.0),Gaussian(m=2.0,V=12.0)]),id=:m1_min_1)
# PriorNode(GaussianDistribution(m=2.0,V=12.0),id=:m2_min_1)
PriorNode(PartitionedDistribution([Gamma(a=1,b=0.01),Gamma(a=1, b=0.01)]),id=:w1_min_1)
# PriorNode(GammaDistribution(a=1,b=0.01),id=:w2_min_1)
PriorNode(Beta(a=2.0,b=2.0),id=:pi_min_1)
TerminalNode(vague(PartitionedDistribution{Gaussian, 2}),id=:m1_n)
# TerminalNode(vague(GaussianDistribution),id=:m2_n)
TerminalNode(vague(PartitionedDistribution{Gamma,2}),id=:w1_n)
# TerminalNode(vague(GammaDistribution),id=:w2_n)
TerminalNode(vague(Beta),id=:pi_n)


Edge(n(:pi_eq).i[1],n(:gm).i[:pi],id=:pi)
Edge(n(:m1_eq).i[1],n(:gm).i[:m],id=:m1)
Edge(n(:w1_eq).i[1],n(:gm).i[:w],id=:w1)
#Edge(n(:m2_eq).i[1],n(:gm).i[:m2],id=:m2)
# Edge(n(:w2_eq).i[1],n(:gm).i[:w2],id=:w2)
Edge(n(:z).i[:out],n(:gm).i[:z],id=:z)
Edge(n(:y).i[:out],n(:gm).i[:x],id=:y)
 Edge(n(:m1_min_1).i[:out], n(:m1_eq).i[:2])
# Edge(n(:m2_min_1).i[:out], n(:m2_eq).i[:2])
 Edge(n(:w1_min_1).i[:out], n(:w1_eq).i[:2])
# Edge(n(:w2_min_1).i[:out], n(:w2_eq).i[:2])
 Edge(n(:pi_min_1).i[:out], n(:pi_eq).i[:2])
 Edge(n(:m1_n).i[:out], n(:m1_eq).i[:3])
# Edge(n(:m2_n).i[:out], n(:m2_eq).i[:3])
 Edge(n(:w1_n).i[:out], n(:w1_eq).i[:3])
# Edge(n(:w2_n).i[:out], n(:w2_eq).i[:3])
 Edge(n(:pi_n).i[:out], n(:pi_eq).i[:3])



 Wrap(n(:pi_n),n(:pi_min_1))
 Wrap(n(:m1_n),n(:m1_min_1))
 #Wrap(n(:m2_n),n(:m2_min_1))
 Wrap(n(:w1_n),n(:w1_min_1))
# Wrap(n(:w2_n),n(:w2_min_1))
draw()

attachReadBuffer(n(:y), deepcopy(y));

m1_est = attachWriteBuffer(n(:m1_eq).i[3])
w1_est = attachWriteBuffer(n(:w1_eq).i[3])
#m2_est = attachWriteBuffer(n(:m2_eq).i[3])
#w2_est = attachWriteBuffer(n(:w2_eq).i[3])
pi_est = attachWriteBuffer(n(:pi_eq).i[3]);

# Specify the variational algorithm for n_its vmp iterations
algo = VariationalBayes(Dict(   eg(:m1) => PartitionedDistribution{Gaussian,2},
                                eg(:w1) => PartitionedDistribution{Gamma,2},
                                eg(:z)  => Bernoulli,
                                eg(:pi) => Beta,
                                eg(:y)  => Gaussian),
                        n_iterations=n_its)

show(algo)
