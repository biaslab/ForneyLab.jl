using ForneyLab
using Distributions

# Initial settings
N                 = [20;20;20]                                # Number of observed samples first clusters
n_its             = 50                                        # Number of vmp iterations
true_mean1        = [10.0,3.0,1.0]                            # Mean of the first cluster
true_variance1    = [3.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]   # Variance of the first cluster
true_mean2        = [4.0,10.0, 2.0]                          # Mean of the second cluster
true_variance2    = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]   # Variance of the second cluster
true_mean3        = [1.0,3.0, 10.0]                          # Mean of the second cluster
true_variance3    = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]   # Variance of the second cluster
d1                = MvNormal(true_mean1, true_variance1)      # Construct the distribution of first cluster
d2                = MvNormal(true_mean2, true_variance2)      # Construct distribution of the second cluster
d3                = MvNormal(true_mean3, true_variance3)      # Construct distribution of the second cluster
y_observations1   = rand(d1,N[1])'                              # Take samples from the first cluster
y_observations2   = rand(d2,N[2])'                              # Take samples from the second cluster
y_observations3   = rand(d3,N[3])'                              # Take samples from the second cluster

permutations      = shuffle(collect(1:sum(N)))
y                 = [y_observations1; y_observations2; y_observations3][permutations,:]        # Mix the samples from the clusters

#Overwrite the Wishart prior with variables that are equal or bigger than 1
function ForneyLab.vague!{dims}(dist::ForneyLab.Wishart{dims})
    dist.V = eye(dims)*4.
    dist.nu = 4.
   return dist
end

ForneyLab.vague{dims}(::Type{ForneyLab.Wishart{dims}}) = ForneyLab.Wishart(V=diageye(dims)*4., nu=4.)

# Build graph
for k=1:(sum(N))
    GaussianMixtureNodePar(id=:gm*k) # s() for symbol concatenation
    EqualityNode(id=:m1_eq*k)
    EqualityNode(id=:w1_eq*k)
    EqualityNode(id=:pi_eq*k)
    TerminalNode(MvDelta(reshape(y[k,:],size(y,2))), id=:y*k) # Observed y values are stored in terminal node values
    PriorNode(ForneyLab.Categorical{length(N)}(1/length(N)*ones(length(N))),id=:z*k)
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

PriorNode(PartitionedDistribution([MvGaussian(m=[5.,1.0, 1.0],V=[10. 0.0 0.0;0.0 10. 0.0; 0.0 0.0 10.0]),MvGaussian(m=[1.0, 5.0, 1.0],V=[10. 0.0 0.0; 0.0 10. 0.0; 0.0 0.0 10.0]), MvGaussian(m=[1.0, 1.0, 5.0],V=[10. 0.0 0.0; 0.0 10. 0.0; 0.0 0.0 10.0])]),id=:m1_start)
PriorNode(PartitionedDistribution([ForneyLab.Wishart(nu=3., V=eye(3)*3.),ForneyLab.Wishart(nu=3., V=eye(3)*3.), ForneyLab.Wishart(nu=3., V=eye(3)*3.)]),id=:w1_start)
PriorNode(ForneyLab.Dirichlet([200.0;200.0;200.0]),id=:pi_start)

Edge(n(:m1_eq*1).i[3],n(:m1_start).i[:out])
Edge(n(:w1_eq*1).i[3],n(:w1_start).i[:out])
Edge(n(:pi_eq*1).i[3],n(:pi_start).i[:out])

TerminalNode(vague(PartitionedDistribution{MvGaussian{length(true_mean1)},length(N)}),id=:m1_end)
TerminalNode(vague(PartitionedDistribution{ForneyLab.Wishart{length(true_mean1)},length(N)}),id=:w1_end)
TerminalNode(vague(ForneyLab.Dirichlet{length(N)}),id=:pi_end)

Edge(n(:m1_eq*(sum(N))).i[2], n(:m1_end))
Edge(n(:w1_eq*(sum(N))).i[2], n(:w1_end))
Edge(n(:pi_eq*(sum(N))).i[2], n(:pi_end))

#Attach write buffers
m1_est = attachWriteBuffer(n(:m1_end).i[:out].partner)
w1_est = attachWriteBuffer(n(:w1_end).i[:out].partner)
pi_est = attachWriteBuffer(n(:pi_end).i[:out].partner)
z_est1=attachWriteBuffer(n(:z*1).i[:out].partner)
z_est2=attachWriteBuffer(n(:z*2).i[:out].partner)
z_est3=attachWriteBuffer(n(:z*3).i[:out].partner)
z_est4=attachWriteBuffer(n(:z*4).i[:out].partner)

# Specify the variational algorithm for n_its vmp iterations
algo = VariationalBayes(Dict(   eg(:m1_e*(1:sum(N))) => PartitionedDistribution{MvGaussian{length(true_mean1)},length(N)},
                                eg(:w1_e*(1:sum(N))) => PartitionedDistribution{ForneyLab.Wishart{length(true_mean1)},length(N)},
                                eg(:z_e*(1:sum(N)))  => ForneyLab.Categorical{length(N)},
                                eg(:pi_e*(1:sum(N))) => ForneyLab.Dirichlet{length(N)},
                                eg(:y_e*(1:sum(N)))  => MvGaussian{length(true_mean1)}),
                        n_iterations=n_its)

show(algo)

run(algo)

#print the true parameters
println("True mean 1: $(true_mean1)")
println("True precision 1: $(inv(true_variance1))")
println("True mean 2: $(true_mean2)")
println("True precision 2: $(inv(true_variance2))")
println("True mean 3: $(true_mean3)")
println("True precision 3: $(inv(true_variance3))")
println("Number of samples: $(sum(N))")

#print the estimated parameters
ensureParameters!(m1_est[end].factors[1], (:m, :V))
ensureParameters!(m1_est[end].factors[2], (:m, :V))
ensureParameters!(m1_est[end].factors[3], (:m, :V))
println("\n----- Online estimation after $(n_its) VMP updates per sample -----")
println("Mean estimate1: $(round(m1_est[end].factors[1].m,2)), with variance $(round((m1_est[end].factors[1].V),2))")
println("Precision estimate1: $(round(w1_est[end].factors[1].nu*w1_est[end].factors[1].V,2))")
println("Mean estimate2: $(round(m1_est[end].factors[2].m,2)), with variance $(round(m1_est[end].factors[2].V,2))")
println("Precision estimate2: $(round(w1_est[end].factors[2].nu*w1_est[end].factors[2].V,2))")
println("Mean estimate3: $(round(m1_est[end].factors[3].m,2)), with variance $(round(m1_est[end].factors[3].V,2))")
println("Precision estimate3: $(round(w1_est[end].factors[3].nu*w1_est[end].factors[3].V,2))")
println("$(z_est1[end])")
println("$(z_est2[end])")
println("$(z_est3[end])")
println("$(z_est4[end])")

#println("ratio of pi $(pi_est[end].a/(pi_est[end].a+pi_est[end].b))")
