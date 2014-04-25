export
    GaussianMessage,
    GeneralMessage

export
    isWellDefined,
    isConsistent,
    ensureMVParametrization!,
    ensureMWParametrization!,
    ensureXiVParametrization!,
    ensureXiWParametrization!

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
    for (key, val) in args
        setfield(self, key, deepcopy(val))
    end

    # In the case of single value V and W, cast V and W to matrix
    self.W = ensureMatrix(self.W)
    self.V = ensureMatrix(self.V)

    # Check parameterizations
    if is(self.m, nothing) && is(self.xi, nothing)
        error("Cannot create GaussianMessage: you should define m or xi or both.")
    end
    if is(self.V, nothing) && is(self.W, nothing)
        error("Cannot create GaussianMessage: you should define V or W or both.")
    end

    return self
end
GaussianMessage() = GaussianMessage(m=[0.0], V=[1.0])

# Methods to check and convert different parametrizations
function isWellDefined(msg::GaussianMessage)
    # Check if msg is not underdetermined
    return !((is(msg.m, nothing) && is(msg.xi, nothing)) ||
             (is(msg.V, nothing) && is(msg.W, nothing)))
end
function isConsistent(msg::GaussianMessage)
    # Check if msg is consistent in case it is overdetermined
    if !is(msg.V, nothing) && !is(msg.W, nothing)
        if maximum(abs(inv(msg.V) - msg.W)) > epsilon
            return false # V and W are not consistent
        end
    end
    if !is(msg.m, nothing) && !is(msg.xi, nothing)
        if !is(msg.V, nothing)
            if maximum(abs(msg.V * msg.xi - msg.m)) > epsilon
                return false # m and xi are not consistent
            end
        else
            if maximum(abs(msg.W * msg.m - msg.xi)) > epsilon
                return false # m and xi are not consistent
            end
        end
    end
    return true # all validations passed
end
function ensureMDefined!(msg::GaussianMessage)
    # Ensure that msg.m is defined, calculate it if needed.
    # An underdetermined msg will throw an exception, we assume msg is well defined.
    msg.m = is(msg.m, nothing) ? ensureVDefined!(msg).V * msg.xi : msg.m
    return msg
end
function ensureXiDefined!(msg::GaussianMessage)
    # Ensure that msg.xi is defined, calculate it if needed.
    # An underdetermined msg will throw an exception, we assume msg is well defined.
    msg.xi = is(msg.xi, nothing) ? ensureWDefined!(msg).W * msg.m : msg.xi
    return msg
end
function ensureVDefined!(msg::GaussianMessage)
    # Ensure that msg.V is defined, calculate it if needed.
    # An underdetermined msg will throw an exception, we assume msg is well defined.
    msg.V = is(msg.V, nothing) ? inv(msg.W) : msg.V
    return msg
end
function ensureWDefined!(msg::GaussianMessage)
    # Ensure that msg.W is defined, calculate it if needed.
    # An underdetermined msg will throw an exception, we assume msg is well defined.
    msg.W = is(msg.W, nothing) ? inv(msg.V) : msg.W
    return msg
end
ensureMVParametrization!(msg::GaussianMessage) = ensureVDefined!(ensureMDefined!(msg))
ensureMWParametrization!(msg::GaussianMessage) = ensureWDefined!(ensureMDefined!(msg))
ensureXiVParametrization!(msg::GaussianMessage) = ensureVDefined!(ensureXiDefined!(msg))
ensureXiWParametrization!(msg::GaussianMessage) = ensureWDefined!(ensureXiDefined!(msg))

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