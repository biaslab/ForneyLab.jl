using ForneyLab

# Initial settings
N               = [100; 100; 100]                               # Number of observed samples first cluster
n_its           = 50                                         # Number of vmp iterations
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
for k=1:(sum(N))
    GaussianMixtureNodePar(id=:gm*k) # s() for symbol concatenation
    EqualityNode(id=:m1_eq*k)
    EqualityNode(id=:w1_eq*k)
    EqualityNode(id=:pi_eq*k)
    TerminalNode(Gaussian(m=y[k],V=10.), id=:y*k) # Observed y values are stored in terminal node values
    PriorNode(Categorical{length(N)}(1/length(N)*ones(length(N))),id=:z*k)
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


PriorNode(PartitionedDistribution([Gaussian(m=25.0,V=12.0),Gaussian(m=10.0,V=12.0),Gaussian(m=120.,V=12.)]),id=:m1_start)
PriorNode(PartitionedDistribution([Gamma(a=1,b=1.),Gamma(a=1, b=1.), Gamma(a=1,b=1.)]),id=:w1_start)
PriorNode(Dirichlet(2.*ones(length(N))),id=:pi_start)

Edge(n(:m1_eq*1).i[3],n(:m1_start).i[:out])
Edge(n(:w1_eq*1).i[3],n(:w1_start).i[:out])
Edge(n(:pi_eq*1).i[3],n(:pi_start).i[:out])

TerminalNode(vague(PartitionedDistribution{Gaussian, length(N)}),id=:m1_end)
TerminalNode(vague(PartitionedDistribution{Gamma,length(N)}),id=:w1_end)
TerminalNode(vague(Dirichlet{length(N)}),id=:pi_end)

Edge(n(:m1_eq*(sum(N))).i[2], n(:m1_end))
Edge(n(:w1_eq*(sum(N))).i[2], n(:w1_end))
Edge(n(:pi_eq*(sum(N))).i[2], n(:pi_end))

# attach write buffers to the wanted variables
m1_est = attachWriteBuffer(n(:m1_end).i[:out].partner)
w1_est = attachWriteBuffer(n(:w1_end).i[:out].partner)
pi_est = attachWriteBuffer(n(:pi_end).i[:out].partner);
x_est1  =  attachWriteBuffer(n(:y1).i[:out].partner);
# Specify the variational algorithm for n_its vmp iterations
msg_types = Dict{Interface,DataType}([n(:gm*i).i[:x] => Gaussian for i=1:(sum(N))])
algo = VariationalBayes(Dict(   eg(:m1_e*(1:(sum(N)))) => PartitionedDistribution{Gaussian,length(N)},
                                eg(:w1_e*(1:(sum(N)))) => PartitionedDistribution{Gamma,length(N)},
                                eg(:z_e*(1:(sum(N))))  => Categorical{length(N)},
                                eg(:pi_e*(1:(sum(N)))) => Dirichlet{length(N)},
                                eg(:y_e*(1:(sum(N))))  => Gaussian),
                        n_iterations=n_its,
                        message_types=msg_types)

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


println(x_est1)
println(y[1])
