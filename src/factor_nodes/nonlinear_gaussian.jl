export NonlinearGaussian

"""
Description:

    NonlinearGaussian node modeling a nonlinear relation with additive Gaussian noise.
    Updates for the NonlinearGaussian node are computed through the unscented transform.
    
    For more details see "On Approximate NonlinearGaussian Gaussian Message Passing on
    Factor Graphs", Petersen et al. 2018.

    f(out, in1, v) = N(out | g(in1), v)

Interfaces:

    1. out
    2. in1
    3. v (covariance)

Construction:

    NonlinearGaussian(out, in1, v, g, id=:my_node)
"""
mutable struct NonlinearGaussian <: SoftFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    g_inv::Union{Function, Nothing} # Inverse of g (optional)
    alpha::Union{Float64, Nothing} # Spread parameter for unscented transform
    dims::Tuple # Dimension of breaker message on input interface

    function NonlinearGaussian(out, in1, v, g::Function; g_inv=nothing, alpha=nothing, dims=(), id=ForneyLab.generateId(NonlinearGaussian))
        @ensureVariables(out, in1, v)
        self = new(id, Vector{Interface}(undef, 3), Dict{Symbol,Interface}(), g, g_inv, alpha, dims)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:v] = self.interfaces[3] = associate!(Interface(self), v)

        return self
    end
end

slug(::Type{NonlinearGaussian}) = "g"

function averageEnergy(::Type{NonlinearGaussian}, 
                       marg_out_in1::ProbabilityDistribution{Multivariate, F},
                       marg_var::ProbabilityDistribution{Univariate},
                       g::Function) where F<:Gaussian

    # Cubature through unscented approximation
    # See Sarkka (2013), "Bayesian Filtering and Smoothing"
    (sigma_points, weights, _) = ForneyLab.sigmaPointsAndWeights(marg_out_in1; alpha=1.0, beta=0.0, kappa=1.0)
    sigma_out = [sigma_point[1] for sigma_point in sigma_points]
    sigma_in1 = [sigma_point[2] for sigma_point in sigma_points]
    t = (sigma_out .- g.(sigma_in1)).^2

    0.5*log(2*pi) +
    0.5*unsafeLogMean(marg_var) +
    0.5*unsafeInverseMean(marg_var)*sum(weights.*t)
end


#--------------------------
# Custom inbounds collector
#--------------------------

function collectAverageEnergyInbounds(node::NonlinearGaussian)
    inbounds = Any[]
    local_posterior_factor_to_region = localPosteriorFactorToRegion(node)

    encountered_posterior_factors = Union{PosteriorFactor, Edge}[] # Keep track of encountered posterior factors
    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        current_posterior_factor = posteriorFactor(node_interface.edge)

        if (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(copy(inbound_interface.node), ProbabilityDistribution)) # Copy Clamp before assembly to prevent overwriting dist_or_msg field
        elseif !(current_posterior_factor in encountered_posterior_factors)
            # Collect marginal entry from marginal dictionary (if marginal entry is not already accepted)
            target = local_posterior_factor_to_region[current_posterior_factor]
            push!(inbounds, current_inference_algorithm.target_to_marginal_entry[target])
        end

        push!(encountered_posterior_factors, current_posterior_factor)
    end

    # Push function to calling signature
    # These functions needs to be defined in the scope of the user
    push!(inbounds, Dict{Symbol, Any}(:g => node.g,
                                      :keyword => false))
    
    return inbounds
end