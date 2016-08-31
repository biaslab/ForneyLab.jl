export
    MvGaussian,
    ensureParameters!,
    isWellDefined,
    isConsistent

"""
Description:

    Encodes a multivariate Gaussian distribution.

Parameters:

    m (“mean”, real vector), V (“variance”, real matrix),
    W (“precision”, real matrix), xi (“weighted mean”, real vector)

Construction:

    MvGaussian(m=[1.0,3.0], V=[2.0, 0.0; 0.0, 2.0])
    MvGaussian(m=[1.0,3.0], V=Diagonal([2.0, 2.0]))

Reference:

    Bishop, 2006; Pattern recognition and machine learning; appendix B
"""
type MvGaussian{dims} <: Multivariate{dims}
    m::Vector{Float64}   # Mean vector
    V::AbstractMatrix{Float64}   # Covariance matrix
    W::AbstractMatrix{Float64}   # Weight matrix
    xi::Vector{Float64}  # Weighted mean vector: xi=W*m

    function MvGaussian(m, V, W, xi)
        (size(m) == size(xi)) || error("Cannot create MvGaussian: m and xi should have the same size")
        (size(V) == size(W)) || error("Cannot create MvGaussian: V and W should have the same size")
        (length(m) == size(V,1) == size(V,2)) || error("Cannot create MvGaussian: inconsistent parameter dimensions")

        if isValid(V)
            (maximum(abs(V)) < realmax(Float64)) || error("Cannot create MvGaussian, covariance matrix V cannot contain Inf.")
            all(abs(diag(V)) .> realmin(Float64)) || error("Cannot create MvGaussian, diagonal of covariance matrix V should be non-zero.")
        end
        if isValid(W)
            (maximum(abs(W)) < realmax(Float64)) || error("Cannot create MvGaussian, precision matrix W cannot contain Inf.")
            all(abs(diag(W)) .> realmin(Float64)) || error("Cannot create MvGaussian, diagonal of precision matrix W should be non-zero.")
        end

        self = new{length(m)}(m, V, W, xi)
        isWellDefined(self) || error("Cannot create MvGaussian, distribution is underdetermined.")

        return self
    end
end

function MvGaussian(; m::Vector{Float64}=[NaN],
                                  V::AbstractMatrix{Float64}=Diagonal([NaN]),
                                  W::AbstractMatrix{Float64}=Diagonal([NaN]),
                                  xi::Vector{Float64}=[NaN])
    # Ensure _m and _xi have the same size
    _m = copy(m)
    _xi = copy(xi)
    if size(_m) != size(_xi)
        if isValid(_m) && !isValid(_xi)
            _xi = fill!(similar(_m), NaN)
        elseif !isValid(_m) && isValid(_xi)
            _m = fill!(similar(_xi), NaN)
        else
            error("m and xi should have the same length")
        end
    end

    # Ensure _V and _W have the same size
    _V = copy(V)
    _W = copy(W)
    if size(_V) != size(_W)
        if isValid(_V) && !isValid(_W)
            _W = fill!(similar(_V), NaN)
        elseif !isValid(_V) && isValid(_W)
            _V = fill!(similar(_W), NaN)
        else
            error("V and W should have the same size")
        end
    end

    return MvGaussian{length(_m)}(_m, _V, _W, _xi)
end

MvGaussian() = MvGaussian(m=[0.0], V=reshape([1.0],1,1))

function pdf(dist::MvGaussian, x::Vector{Float64})
    k = size(dist.V, 1)
    if isValid(dist.W)
        if isValid(dist.m)
            C = (2.0*pi)^(-k/2.0) * det(dist.W)^(0.5)
            d = x - dist.m
            return (C * exp(-0.5 * d' * dist.W * d))[1]
        elseif isValid(dist.xi)
            W = dist.W; ξ = dist.xi
            C = (2.0*pi)^(-k/2.0) * det(W)^(0.5)
            return (C * exp(-0.5*ξ'*pinv(W)*ξ) * exp(-0.5*x'*W*x + x'*ξ))[1]
        else
            error("Cannot evaluate pdf for underdetermined MvGaussian")
        end
    elseif isValid(dist.V)
        ensureParameter!(dist, Val{:m})
        C = (2.0*pi)^(-k/2.0) * det(dist.V)^(-0.5)
        d = x - dist.m
        return (C * exp(-0.5 * d' * pinv(dist.V) * d))[1]
    else
        error("Cannot evaluate pdf for underdetermined MvGaussian")
    end
end

function vague!{dims}(dist::MvGaussian{dims})
    dist.m = zeros(dims)
    dist.V = huge*diageye(dims)
    invalidate!(dist.W)
    invalidate!(dist.xi)
    return dist
end

vague{dims}(::Type{MvGaussian{dims}}) = MvGaussian(m=zeros(dims), V=huge*diageye(dims))

function format(dist::MvGaussian)
    if isValid(dist.m) && isValid(dist.V)
        return "N(m=$(format(dist.m)), V=$(format(dist.V)))"
    elseif isValid(dist.m) && isValid(dist.W)
        return "N(m=$(format(dist.m)), W=$(format(dist.W)))"
    elseif isValid(dist.xi) && isValid(dist.W)
        return "N(ξ=$(format(dist.xi)), W=$(format(dist.W)))"
    elseif isValid(dist.xi) && isValid(dist.V)
        return "N(ξ=$(format(dist.xi)), V=$(format(dist.V)))"
    else
        return "N(underdetermined)"
    end
end

show(io::IO, dist::MvGaussian) = println(io, format(dist))

function prod!{dims}(x::MvGaussian{dims}, y::MvGaussian{dims}, z::MvGaussian{dims}=MvGaussian(m=zeros(dims), V=eye(dims)))
    # Multiplication of 2 multivariate Gaussian PDFs: p(z) = p(x) * p(y)
    # Equations from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
    if isValid(x.xi) && isValid(x.W) && isValid(y.xi) && isValid(y.W)
        # Use (xi,W)x(xi,W) parametrization
        invalidate!(z.m)
        invalidate!(z.V)
        z.W  = x.W + y.W
        z.xi = x.xi + y.xi
    elseif isValid(x.m) && isValid(x.W) && isValid(y.m) && isValid(y.W)
        # Use (m,W)x(m,W) parametrization
        z.m  = cholinv(x.W+y.W) * (x.W*x.m + y.W*y.m)
        invalidate!(z.V)
        z.W  = x.W + y.W
        invalidate!(z.xi)
    elseif isValid(x.xi) && isValid(x.V) && isValid(y.xi) && isValid(y.V)
        # Use (xi,V)x(xi,V) parametrization
        invalidate!(z.m)
        z.V  = x.V * cholinv(x.V+y.V) * y.V
        invalidate!(z.W)
        z.xi = x.xi + y.xi
    else
        # Convert to (xi,W)x(xi,W) parametrization
        ensureParameters!(x, (:xi, :W))
        ensureParameters!(y, (:xi, :W))
        invalidate!(z.m)
        invalidate!(z.V)
        z.W  = x.W + y.W
        z.xi = x.xi + y.xi
    end

    return z
end

@symmetrical function prod!{dims}(x::MvGaussian{dims}, y::MvDelta{Float64,dims}, z::MvDelta{Float64,dims}=MvDelta(y.m))
    # Product of multivariate Gaussian PDF and MvDelta
    (z.m == y.m) || (z.m[:] = y.m)

    return z
end

@symmetrical function prod!{dims}(x::MvGaussian{dims}, y::MvDelta{Float64,dims}, z::MvGaussian{dims})
    # Product of multivariate Gaussian PDF and MvDelta, force result to be MvGaussian
    z.m[:] = y.m
    z.V = tiny.*eye(dims)
    invalidate!(z.xi)
    invalidate!(z.W)

    return z
end

@symmetrical function prod!{dims}(::Void, y::MvDelta{Float64,dims}, z::MvGaussian{dims})
    # Product of an unknown with MvDelta, force result to be MvGaussian
    z.m[:] = y.m
    z.V = tiny.*eye(dims)
    invalidate!(z.xi)
    invalidate!(z.W)

    return z
end

unsafeMean(dist::MvGaussian) = ensureParameter!(dist, Val{:m}).m

unsafeVar(dist::MvGaussian) = diag(ensureParameter!(dist, Val{:V}).V)

unsafeCov(dist::MvGaussian) = ensureParameter!(dist, Val{:V}).V

function isProper(dist::MvGaussian)
    if isWellDefined(dist)
        param = isValid(dist.V) ? dist.V : dist.W
        if isRoundedPosDef(param)
            return true
        end
    end

    return false
end

function sample(dist::MvGaussian)
    isProper(dist) || error("Cannot sample from improper distribution")
    ensureParameters!(dist, (:m, :V))
    return (dist.V^0.5)*randn(length(dist.m)) + dist.m
end

# Methods to check and convert different parametrizations
function isWellDefined(dist::MvGaussian)
    # Check if dist is not underdetermined
    if !((isValid(dist.m) || isValid(dist.xi)) &&
         (isValid(dist.V) || isValid(dist.W)))
        return false
    end
    dimensions=0
    for field in [:m, :xi, :V, :W]
        if isValid(getfield(dist, field))
            if dimensions>0
                if maximum(size(getfield(dist, field))) != dimensions
                    return false
                end
            else
                dimensions = size(getfield(dist, field), 1)
            end
        end
    end
    return true
end

function isConsistent(dist::MvGaussian)
    # Check if dist is consistent in case it is overdetermined
    if isValid(dist.V) && isValid(dist.W)
        V_W_consistent = false
        try
           V_W_consistent = isApproxEqual(cholinv(dist.V), dist.W)
        catch
            try
                V_W_consistent = isApproxEqual(cholinv(dist.W), dist.V)
            catch
                error("Cannot check consistency of MvGaussian because both V and W are non-invertible.")
            end
        end
        if !V_W_consistent
            return false # V and W are not consistent
        end
    end
    if isValid(dist.m) && isValid(dist.xi)
        if isValid(dist.V)
            if isApproxEqual(dist.V * dist.xi, dist.m) == false
                return false # m and xi are not consistent
            end
        else
            if isApproxEqual(dist.W * dist.m, dist.xi) == false
                return false # m and xi are not consistent
            end
        end
    end
    return true # all validations passed
end

# In all ensureParameter! methods we check if the required parameter defined and, if not, calculate it.
# We assume that the distribution is well-defined, otherwise we would've gotten the message upon creating.

function ensureParameters!(dist::MvGaussian, params::Tuple{Symbol, Vararg{Symbol}})
    for param in params
        ensureParameter!(dist, Val{param})
    end
    return dist
end

function ensureParameter!(dist::MvGaussian, param::Type{Val{:m}})
    if !isValid(dist.m)
        dist.m = ensureParameter!(dist, Val{:V}).V * dist.xi
    end
    return dist
end

function ensureParameter!(dist::MvGaussian, param::Type{Val{:xi}})
    if !isValid(dist.xi)
        dist.xi = ensureParameter!(dist, Val{:W}).W * dist.m
    end
    return dist
end

function ensureParameter!(dist::MvGaussian, param::Type{Val{:W}})
    if !isValid(dist.W)
        dist.W = cholinv(dist.V)
    end
    return dist
end

function ensureParameter!(dist::MvGaussian, param::Type{Val{:V}})
    if !isValid(dist.V)
        dist.V = cholinv(dist.W)
    end
    return dist
end

function ==(x::MvGaussian, y::MvGaussian)
    if is(x, y)
        return true
    end
    if !isWellDefined(x) || !isWellDefined(y)
        return false
    end
    # Check m or xi
    if isValid(x.m) && isValid(y.m)
        (length(x.m)==length(y.m)) || return false
        isApproxEqual(x.m,y.m) || return false
    elseif isValid(x.xi) && isValid(y.xi)
        (length(x.xi)==length(y.xi)) || return false
        isApproxEqual(x.xi,y.xi) || return false
    else
        ensureParameter!(x, Val{:m}); ensureParameter!(y, Val{:m});
        (length(x.m)==length(y.m)) || return false
        isApproxEqual(x.m,y.m) || return false
    end

    # Check V or W
    if isValid(x.V) && isValid(y.V)
        (length(x.V)==length(y.V)) || return false
        isApproxEqual(x.V,y.V) || return false
    elseif isValid(x.W) && isValid(y.W)
        (length(x.W)==length(y.W)) || return false
        isApproxEqual(x.W,y.W) || return false
    else
        ensureParameter!(x, Val{:V}); ensureParameter!(y, Val{:V});
        (length(x.V)==length(y.V)) || return false
        isApproxEqual(x.V,y.V) || return false
    end

    return true
end

# Convert Delta -> MvGaussian
# NOTE: this introduces a small error because the variance is set >0
convert(::Type{MvGaussian}, delta::MvDelta{Float64}) = MvGaussian(m=delta.m, V=tiny*diageye(length(delta.m)))
convert{TD<:MvDelta{Float64}, TG<:MvGaussian}(::Type{Message{TG}}, msg::Message{TD}) = Message(MvGaussian(m=msg.payload.m, V=tiny*diageye(length(msg.payload.m))))

# Convert Gaussian -> MvGaussian
convert(::Type{MvGaussian}, d::Gaussian) = MvGaussian(m=[d.m], V=d.V*eye(1), W=d.W*eye(1), xi=[d.xi])

# Convert MvGaussian -> Gaussian
function convert(::Type{Gaussian}, d::MvGaussian{1})
    Gaussian(m=d.m[1], V=d.V[1,1], W=d.W[1,1], xi=d.xi[1])
end

# Entropy functional
function H{dims}(dist::MvGaussian{dims})
    ensureParameters!(dist, (:m, :V))

    return  0.5*log(det(dist.V)) +
            (dims/2)*log(2*pi) +
            (dims/2)    
end