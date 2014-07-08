export # types
    GammaMessage,
    InverseGammaMessage,
    GaussianMessage,
    GeneralMessage

export # functions
    ensureMVParametrization!,
    ensureMWParametrization!,
    ensureXiVParametrization!,
    ensureXiWParametrization!,
    isWellDefined,
    isConsistent,
    uninformative,
    getOrAssign

function getOrAssign(interface::Interface, assign_type::DataType)
    # Looks for a message on interface.
    # When no message is present, it sets and returns a standard message.
    # Otherwise it returns the present message.

    if interface.message == nothing
        msg_out = assign_type() # Initialize a message if there is none
        interface.message = msg_out
    else
        msg_out = interface.message
    end
    return msg_out
end

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
    m::Union(Array{Float64, 1}, Nothing)    # Mean vector
    V::Union(Array{Float64}, Nothing)       # Covariance matrix
    W::Union(Array{Float64}, Nothing)       # Weight matrix
    xi::Union(Array{Float64, 1}, Nothing)   # Weighted mean vector: xi=W*m
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
    if !isWellDefined(self)
        error("Cannot create GaussianMessage, parameterization is underdetermined.")
    end

    return self
end
GaussianMessage() = GaussianMessage(m=[0.0], V=[1.0])

uninformative(msg_type::Type{GaussianMessage}) = GaussianMessage(m=[0.0], V=[1000.0])

function show(io::IO, msg::GaussianMessage)
    println(io, "GaussianMessage")
    print(io, "m  = ")
    show(io, msg.m)
    print(io, "\nV  = ")
    show(io, msg.V)
    print(io, "\nW  = ")
    show(io, msg.W)
    print(io, "\nxi = ")
    show(io, msg.xi)
    print(io, "\n")
end

Base.mean(msg::GaussianMessage) = ensureMDefined!(msg).m
Base.var(msg::GaussianMessage) = diag(ensureVDefined!(msg).V, 0)

# Methods to check and convert different parametrizations
function isWellDefined(msg::GaussianMessage)
    # Check if msg is not underdetermined
    if ((is(msg.m, nothing) && is(msg.xi, nothing)) ||
        (is(msg.V, nothing) && is(msg.W, nothing)))
        return false
    end
    dimensions=0
    for field in [:m, :xi, :V, :W]
        if !is(getfield(msg, field), nothing)
            if dimensions>0
                if maximum(size(getfield(msg, field)))!=dimensions
                    return false
                end
            else
                dimensions = maximum(size(getfield(msg, field)))
            end
        end
    end
    return true
end
function isConsistent(msg::GaussianMessage)
    # Check if msg is consistent in case it is overdetermined
    if !is(msg.V, nothing) && !is(msg.W, nothing)
        V_W_consistent = false
        try
           V_W_consistent = isApproxEqual(inv(msg.V), msg.W)
        catch
            try
                V_W_consistent = isApproxEqual(inv(msg.W), msg.V)
            catch
                error("Cannot check consistency of GaussianMessage because both V and W are non-invertible.")
            end
        end
        if !V_W_consistent
            return false # V and W are not consistent
        end
    end
    if !is(msg.m, nothing) && !is(msg.xi, nothing)
        if !is(msg.V, nothing)
            if isApproxEqual(msg.V * msg.xi, msg.m) == false
                return false # m and xi are not consistent
            end
        else
            if isApproxEqual(msg.W * msg.m, msg.xi) == false
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
    try
        msg.V = is(msg.V, nothing) ? inv(msg.W) : msg.V
    catch
        error("Cannot calculate V of GaussianMessage because W is not invertible.")
    end
    return msg
end
function ensureWDefined!(msg::GaussianMessage)
    # Ensure that msg.W is defined, calculate it if needed.
    # An underdetermined msg will throw an exception, we assume msg is well defined.
    try
        msg.W = is(msg.W, nothing) ? inv(msg.V) : msg.W
    catch
        error("Cannot calculate W of GaussianMessage because V is not invertible.")
    end
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
        self = new()
        self.value = deepcopy(value) # Make a copy instead of referencing
        return self
    end
end
GeneralMessage() = GeneralMessage(1.0)

uninformative(msg_type::Type{GeneralMessage}) = GeneralMessage(1.0)

function show(io::IO, msg::GeneralMessage)
    print(io, "GeneralMessage with value = ")
    show(io, msg.value)
    print("\n")
end

############################################
# GammaMessage
############################################
# Description:
#   Encodes a gamma PDF.
#   Pamameters: scalars a (shape) and b (rate).
############################################
type GammaMessage <: Message
    a::Float64 # shape
    b::Float64 # rate
    GammaMessage(; a=1.0, b=1.0) = new(a, b)
end

uninformative(msg_type::Type{GammaMessage}) = GammaMessage(a=0.999, b=0.001)
Base.mean(msg::GammaMessage) = msg.a / msg.b
Base.var(msg::GammaMessage) = msg.a / (msg.b^2)

############################################
# InverseGammaMessage
############################################
# Description:
#   Encodes an inverse gamma PDF.
#   Pamameters: scalars a (shape) and b (rate).
############################################
type InverseGammaMessage <: Message
    a::Float64 # shape
    b::Float64 # rate
    InverseGammaMessage(; a=1.0, b=1.0) = new(a, b)
end

uninformative(msg_type::Type{InverseGammaMessage}) = InverseGammaMessage(a=-0.999, b=0.001)
function Base.mean(msg::InverseGammaMessage)
    if msg.a > 1.0
        return msg.b / (msg.a - 1)
    else
        return NaN
    end
end
function Base.var(msg::InverseGammaMessage)
    if msg.a > 2.0
        return (msg.b^2) / ( ( (msg.a-1)^2 ) * (msg.a-2) )
    else
        return NaN
    end
end
function show(io::IO, msg::Union(GammaMessage, InverseGammaMessage))
    println(io, typeof(msg))
    println(io, "a = $(msg.a) (shape)")
    println(io, "b = $(msg.b) (rate)")
end