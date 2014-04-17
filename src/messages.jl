# TODO: Investigate if it's better to define messages as immutable.

export
    GaussianMessage,
    GeneralMessage

############################################
# GaussianMessage
############################################
# Description:
#   Encodes a Gaussian PDF.
#   Define (mean (m) or weighted mean (xi))
#   and (covariance (V) or precision (W)).
#   These result in the same PDF:
#    msg = GaussianMessage(m=[1.0], V=[2.0])
#    msg = GaussianMessage(m=[1.0], W=[0.5])
#    msg = GaussianMessage(xi=[0.5], W=[0.5])
#    msg = GaussianMessage(xi=[0.5], V=[2.0])
#   m and xi are 1d arrays, V and W are 2d
############################################
type GaussianMessage <: Message
    m::Union(Array{Float64, 1}, Nothing)     # Mean vector
    V::Union(Array{Float64}, Nothing)     # Covariance matrix
    W::Union(Array{Float64}, Nothing)     # Weight matrix
    xi::Union(Array{Float64, 1}, Nothing)    # Weighted mean vector: xi=W*m
end
function GaussianMessage(;args...)
    self = GaussianMessage(nothing, nothing, nothing, nothing)
    for (key,val) in args
        setfield(self, key, deepcopy(val))
    end
    if is(self.m, nothing) && is(self.xi, nothing)
        error("Cannot create GaussianMessage: you should define m or xi or both.")
    end
    if is(self.V, nothing) && is(self.W, nothing)
        error("Cannot create GaussianMessage: you should define V or W or both.")
    end
    if typeof(self.xi)==Array && is(self.W, nothing)
        error("Cannot create GaussianMessage: you should also define W if you use xi")
    end

    return self
end
GaussianMessage() = GaussianMessage(m=[0.0], V=[1.0])

############################################
# GeneralMessage
############################################
# Description:
#   Simply holds an arbitrary object.
#   Useful for example for passing parameters.
############################################
type GeneralMessage <: Message
    value
    function GeneralMessage(value)
        # In case value is mutable, copy it instead of just referencing.
        if (isimmutable(value))
            self = new(value)
        else
            self = new()
            self.value = deepcopy(value)
        end
        return self
    end
end
GeneralMessage() = GeneralMessage(1.0)