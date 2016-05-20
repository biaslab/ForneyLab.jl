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

    function Mixture(components::Vector{dtype}, weights::Vector{Float64})
        (length(components) > 0) || error("Mixture should contain at least one component")
        (all(map(typeof, components) .== typeof(components[1]))) || error("All mixture components should have the same type")

        return new{dtype}(components, weights)
    end
end

Mixture{dtype<:ProbabilityDistribution}(components::Vector{dtype}, weights::Vector{Float64}) = Mixture{typeof(components[1])}(components, weights)

Mixture() = Mixture([Gaussian()], [1.0])

dimensions{dtype}(dist::Mixture{dtype}) = dimensions(dtype)

dimensions{T<:Mixture}(dist_type::Type{T}) = dimensions(dist_type.parameters[1])

function vague!(dist::Mixture)
    for component in dist.components
        vague!(component)
    end
    fill!(dist.weights, 1./length(dist.components))

    return dist
end

vague{T<:Mixture}(dist_type::Type{T}) = Mixture([vague(T.parameters[1])], [1.0])

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

function ==(x::Mixture, y::Mixture)
    is(x,y) && return true # same object
    (typeof(x) == typeof(y)) || return false # different component types
    (length(x.components)==length(y.components)) || return false
    if (x.weights==y.weights) && (x.components==y.components)
        return true
    end

    # Check if x and y are identical except for the ordering of the components
    n = length(x.components)
    permutation_candidates = collect(1:n)
    for i=1:n
        # Find component in y that matches x.components[i]
        matched = false
        for p=1:length(permutation_candidates)
            j = permutation_candidates[p]
            if (x.weights[i]==y.weights[j]) && (x.components[i]==y.components[j])
                matched = true
                deleteat!(permutation_candidates, p)
                break
            end
        end
        matched || return false # x.components[i] could not be matched to a component of y
    end

    return true
end

Base.mean(dist::Mixture) = sum(map(mean, dist.components) .* dist.weights)

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