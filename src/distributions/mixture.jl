export Mixture, prune!

"""
Description:

    A weighted mixture of distributions of the same type.

Pamameters:

    components: vector of component distributions.
    weights: vector of component weights. Should sum to 1.

Construction:

    Mixture([Gaussian();Gaussian()], [0.5;0.5])
"""
type Mixture{dtype<:ProbabilityDistribution} <: ProbabilityDistribution
    components::Vector{dtype}
    weights::Vector{Float64}
end

Mixture() = Mixture([Gaussian()], [1.0])

dimensions{dtype}(dist::Mixture{dtype}) = dimensions(dtype)

dimensions{dtype}(dist_type::Type{Mixture{dtype}}) = dimensions(dtype)

function vague!(dist::Mixture)
    for component in dist.components
        vague!(component)
    end
    fill!(dist.weights, 1./length(dist.components))

    return dist
end

function format(dist::Mixture)
    str = "$(typeof(dist))\n"
    str *= "  weights: $(format(dist.weights))\n"
    str *= "  components: $(format(dist.components))"

    return str
end

show(io::IO, dist::Mixture) = println(io, format(dist))

function isProper(dist::Mixture)
    isApproxEqual(sum(dist.weights), 1.0) || return false
    return all(map(isProper, dist.components))
end

"""
Change the number of components in a Mixture
"""
function Base.resize!{dtype}(dist::Mixture{dtype}, num_components::Int64)
    current_length = length(dist.components)
    if num_components < current_length
        # Delete components
        return deleteat!(dist, num_components:current_length)
    elseif num_components > current_length
        # Add components
        dist.weights = [dist.weights; zeros(num_components-current_length)]
        dist.components = [dist.components; [vague(dtype) for i=1:num_components-current_length]]
        return dist
    else
        return dist # nothing to do
    end
end

"""
Delete one or more components from a Mixture
"""
function Base.deleteat!(dist::Mixture, args...)
    deleteat!(dist.components, args...)
    deleteat!(dist.weights, args...)
    norm = sum(dist.weights)
    (norm > 0.) || error("Cannot normalize Mixture weights after deleting components")
    dist.weights /= norm

    return dist
end


"""
`prune(dist::Mixture; max_num_components=0, min_weight=1e-6)`

Reduce the number of components of a Mixture.
"""
function prune!(dist::Mixture; max_num_components::Int64=0, min_weight::Float64=1e-6)
    # Delete components untill all weights are above the threshold
    while length(dist.components) > 1
        (smallest_weight, idx) = findmin(dist.weights)
        (smallest_weight < min_weight) || break
        deleteat!(dist, idx)
    end

    # Delete additional components until max_num_components is satisfied
    # For now, delete the components with the smallest weights.
    # This approach might be improved by picking components such that the
    # KL divergence between the pruned mixture and the original one is minimized.
    if (max_num_components > 0) && (length(dist.components) > max_num_components)
        # Find indexes of components to remove
        components_to_remove = sortperm(dist.weights, rev=true)[max_num_components+1:end]
        deleteat!(dist, components_to_remove)
    end

    return dist
end