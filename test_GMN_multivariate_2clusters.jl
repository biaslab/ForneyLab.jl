using ForneyLab
using Distributions

# Initial settings
N1                = 50                                        # Number of observed samples first clusters
N2                = 100                                        # Number of observed samples second cluster
n_its             = 50                                        # Number of vmp iterations
true_mean1        = [20.0,5.0,6.0]                            # Mean of the first cluster
true_variance1    = [3.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]   # Variance of the first cluster
true_mean2        = [3.0,10.0, 20.0]                          # Mean of the second cluster
true_variance2    = [1.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 1.0]   # Variance of the second cluster
d1                = MvNormal(true_mean1, true_variance1)      # Construct the distribution of first cluster
d2                = MvNormal(true_mean2, true_variance2)      # Construct distribution of the second cluster
y_observations1   = rand(d1,N1)'                                # Take samples from the first cluster
y_observations2   = rand(d2,N2)'                                # Take samples from the second cluster                            # Take samples from the second cluster

permutations      = shuffle(collect(1:N1+N2))
y                 = [y_observations1; y_observations2][permutations,:]         # Mix the samples from the clusters

#Overwrite the Wishart prior with variables that are equal or bigger than 1
function ForneyLab.vague!{dims}(dist::ForneyLab.Wishart{dims})
    dist.V = eye(dims)
    dist.nu = 1.
   return dist
end

ForneyLab.vague{dims}(::Type{ForneyLab.Wishart{dims}}) = ForneyLab.Wishart(V=diageye(dims), nu=1.)


# Build graph
for k=1:(N1+N2)
    GaussianMixtureNodePar(id=:gm*k) # s() for symbol concatenation
    EqualityNode(id=:m1_eq*k)
    EqualityNode(id=:w1_eq*k)
    EqualityNode(id=:pi_eq*k)
    TerminalNode(MvDelta(reshape(y[k,:],size(y,2))), id=:y*k) # Observed y values are stored in terminal node values
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


PriorNode(PartitionedDistribution([MvGaussian(m=[23.,10.0, 1.0],V=[10. 0.0 0.0;0.0 10. 0.0; 0.0 0.0 10.0]),MvGaussian(m=[9.0, 10.0, 15.0],V=[10. 0.0 0.0; 0.0 10. 0.0; 0.0 0.0 10.0])]),id=:m1_start)
PriorNode(PartitionedDistribution([ForneyLab.Wishart(nu=3., V=eye(3)/3.),ForneyLab.Wishart(nu=3., V=eye(3)/3.)]),id=:w1_start)
PriorNode(ForneyLab.Beta(a=2.,b=2.),id=:pi_start)

Edge(n(:m1_eq*1).i[3],n(:m1_start).i[:out])
Edge(n(:w1_eq*1).i[3],n(:w1_start).i[:out])
Edge(n(:pi_eq*1).i[3],n(:pi_start).i[:out])

TerminalNode(vague(PartitionedDistribution{MvGaussian{length(true_mean1)},2}),id=:m1_end)
TerminalNode(vague(PartitionedDistribution{ForneyLab.Wishart{length(true_mean1)},2}),id=:w1_end)
TerminalNode(vague(ForneyLab.Beta),id=:pi_end)

Edge(n(:m1_eq*(N1+N2)).i[2], n(:m1_end))
Edge(n(:w1_eq*(N1+N2)).i[2], n(:w1_end))
Edge(n(:pi_eq*(N1+N2)).i[2], n(:pi_end))


#Attach write buffers
m1_est = attachWriteBuffer(n(:m1_end).i[:out].partner)
w1_est = attachWriteBuffer(n(:w1_end).i[:out].partner)
pi_est = attachWriteBuffer(n(:pi_end).i[:out].partner)

# Specify the variational algorithm for n_its vmp iterations
msg_types = Dict{Interface,DataType}([n(:gm*i).i[:x] => MvGaussian{length(true_mean1)} for i=1:(N1+N2)])

algo = VariationalBayes(Dict(   eg(:m1_e*(1:(N1+N2))) => PartitionedDistribution{MvGaussian{length(true_mean1)},2},
                                eg(:w1_e*(1:(N1+N2))) => PartitionedDistribution{ForneyLab.Wishart{length(true_mean1)},2},
                                eg(:z_e*(1:(N1+N2)))  => ForneyLab.Bernoulli,
                                eg(:pi_e*(1:(N1+N2))) => ForneyLab.Beta,
                                eg(:y_e*(1:(N1+N2)))  => MvGaussian{length(true_mean1)}),
                        n_iterations=n_its,
                        message_types=msg_types)

show(algo)

run(algo);

#print the true parameters
# println("True mean: $(true_mean1)")
# println("True precision: $(inv(true_variance1))")
# println("True mean: $(true_mean2)")
# println("True precision: $(inv(true_variance2))")
# println("Number of samples: $(N1+N2)")

#print the estimated parameters
ensureParameters!(m1_est[end].factors[1], (:m, :V))
ensureParameters!(m1_est[end].factors[2], (:m, :V))
println("\n----- Online estimation after $(n_its) VMP updates per sample -----")
println("Mean estimate1: $(round(m1_est[end].factors[1].m,2)), with variance $(round((m1_est[end].factors[1].V),2))")
println("Precision estimate1: $(round(w1_est[end].factors[1].nu*w1_est[end].factors[1].V,2))")
println("Mean estimate2: $(round(m1_est[end].factors[2].m,2)), with variance $(round(m1_est[end].factors[2].V,2))")
println("Precision estimate2: $(round(w1_est[end].factors[2].nu*w1_est[end].factors[2].V,2))")
println("ratio of pi $(pi_est[end].a/(pi_est[end].a+pi_est[end].b))")
