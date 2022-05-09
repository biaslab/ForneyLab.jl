export Wishart

"""
Description:

    A Wishart node:

    f(out,v,nu) = W(out|v, nu) = B(v, nu) |out|^{(nu - D - 1)/2} exp(-1/2 tr(v^{-1} out))

Interfaces:

    1. out
    2. v (scale matrix)
    3. nu (degrees of freedom)

Construction:

    Wishart(out, v, nu, id=:some_id)
"""
mutable struct Wishart <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Wishart(out, v, nu; id=generateId(Wishart))
        @ensureVariables(out, v, nu)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:v] = self.interfaces[2] = associate!(Interface(self), v)
        self.i[:nu] = self.interfaces[3] = associate!(Interface(self), nu)

        return self
    end
end

slug(::Type{Wishart}) = "W"

format(dist::Distribution{MatrixVariate, Wishart}) = "$(slug(Wishart))(v=$(format(dist.params[:v])), nu=$(format(dist.params[:nu])))"

Distribution(::Type{MatrixVariate}, ::Type{Wishart}; v=mat(1.0), nu=1.0) = Distribution{MatrixVariate, Wishart}(Dict(:v=>v, :nu=>nu))
Distribution(::Type{Wishart}; v=mat(1.0), nu=1.0) = Distribution{MatrixVariate, Wishart}(Dict(:v=>v, :nu=>nu))

dims(dist::Distribution{MatrixVariate, Wishart}) = size(dist.params[:v])

vague(::Type{Wishart}, dims::Tuple{Int64, Int64}) = Distribution(MatrixVariate, Wishart, v=huge*diageye(dims[1]), nu=Float64(dims[1]))

vague(::Type{Union{Gamma, Wishart}}, dims::Tuple{Int64, Int64}) = vague(Wishart, dims)
vague(::Type{Union{Gamma, Wishart}}, dims::Tuple) = vague(Gamma) # Univariate fallback

unsafeMean(dist::Distribution{MatrixVariate, Wishart}) = dist.params[:nu]*dist.params[:v] # unsafe mean

function unsafeDetLogMean(dist::Distribution{MatrixVariate, Wishart})
    d = dims(dist)[1]
    sum([digamma.(0.5*(dist.params[:nu] + 1 - i)) for i = 1:d]) +
    d*log(2) +
    logdet(dist.params[:v])
end

function unsafeVar(dist::Distribution{MatrixVariate, Wishart}) # unsafe variance
    d = dims(dist)[1]
    M = fill!(similar(Matrix(dist.params[:v])), NaN)
    for i = 1:d
        for j = 1:d
            M[i, j] = dist.params[:nu]*(dist.params[:v][i, j]^2 + dist.params[:v][i, i]*dist.params[:v][j, j])
        end
    end
    return M
end

function logPdf(dist::Distribution{MatrixVariate, Wishart},x)
    d = dims(dist)[1]
    0.5*((dist.params[:nu]-d-1)*logdet(x) - tr(inv(dist.params[:v])*x) - dist.params[:nu]*d*log(2) - dist.params[:nu]*logdet(dist.params[:v])) - logmvgamma(d,0.5*dist.params[:nu])
end

function isProper(dist::Distribution{MatrixVariate, Wishart})
    (size(dist.params[:v], 1) == size(dist.params[:v], 2)) || return false
    (dist.params[:nu] > size(dist.params[:v], 1) - 1) || return false
    isRoundedPosDef(dist.params[:v]) || return false
    return true
end

function prod!( x::Distribution{MatrixVariate, Wishart},
                y::Distribution{MatrixVariate, Wishart},
                z::Distribution{MatrixVariate, Wishart}=Distribution(MatrixVariate, Wishart, v=mat(1.0), nu=1.0))

    d = dims(x)[1]
    z.params[:v] = x.params[:v] * cholinv(x.params[:v] + y.params[:v]) * y.params[:v]
    z.params[:nu] = x.params[:nu] + y.params[:nu] - d - 1.0

    return z
end

@symmetrical function prod!(x::Distribution{MatrixVariate, Wishart},
                            y::Distribution{MatrixVariate, PointMass},
                            z::Distribution{MatrixVariate, PointMass}=Distribution(MatrixVariate, PointMass, m=mat(NaN)))

    isRoundedPosDef(y.params[:m]) || error("PointMass location $(y.params[:m]) should be positive definite")
    z.params[:m] = deepcopy(y.params[:m])

    return z
end

function naturalParams(dist::Distribution{MatrixVariate, Wishart})
    d = dims(dist)[1]
    return vcat(-0.5*vec(cholinv(dist.params[:v])), 0.5*(dist.params[:nu]-d-1))
end

function standardDistribution(V::Type{MatrixVariate}, F::Type{Wishart}; η::Vector)
    d = Int(sqrt(length(η) - 1))
    η_1 = reshape(η[1:end-1], d, d)
    η_2 = η[end]
    return Distribution(V, F, v=cholinv(-2.0*η_1), nu=2*η_2+d+1)
end

function logNormalizer(::Type{MatrixVariate}, ::Type{Wishart}; η::Vector)
    d = Int(sqrt(length(η) - 1))
    η_1 = reshape(η[1:end-1], d, d)
    η_2 = η[end]
    return -(η_2+(d+1)/2)*logdet(-η_1) + logmvgamma(d, η_2+(d+1)/2)
end

logPdf(V::Type{MatrixVariate}, F::Type{Wishart}, x::Matrix; η::Vector) = vcat(vec(x), logdet(x))'*η - logNormalizer(V, F; η=η)

# Entropy functional
function differentialEntropy(dist::Distribution{MatrixVariate, Wishart})
    d = dims(dist)[1]
    0.5*(d + 1.0)*logdet(dist.params[:v]) +
    0.5*d*(d + 1.0)*log(2) +
    0.25*d*(d - 1.0)*log(pi) +
    sum([labsgamma(0.5*(dist.params[:nu] + 1.0 - i)) for i=1:d]) -
    0.5*(dist.params[:nu] - d - 1.0) * sum([digamma.(0.5*(dist.params[:nu] + 1.0 - i)) for i=1:d]) +
    0.5*dist.params[:nu]*d
end

# Average energy functional
function averageEnergy(::Type{Wishart}, marg_out::Distribution{MatrixVariate}, marg_v::Distribution{MatrixVariate}, marg_nu::Distribution{Univariate, PointMass})
    d = dims(marg_out)[1]
    0.5*marg_nu.params[:m]*unsafeDetLogMean(marg_v) +
    0.5*marg_nu.params[:m]*d*log(2) +
    0.25*d*(d - 1.0)*log(pi) +
    sum([labsgamma(0.5*(marg_nu.params[:m] + 1.0 - i)) for i=1:d]) -
    0.5*(marg_nu.params[:m] - d - 1.0)*unsafeDetLogMean(marg_out) +
    0.5*tr(unsafeInverseMean(marg_v)*unsafeMean(marg_out))
end
