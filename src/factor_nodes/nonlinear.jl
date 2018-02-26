export Nonlinear

"""
Description:

    Nonlinear node modeling a nonlinear relation with added Gaussian noise.
    Updates for the nonlinear node are computed through local linearization.

    f(out, in1, w) = 𝒩(out| g(in1), w^{-1})

Interfaces:

    1. out
    2. in1
    5. w (precision)

Construction:

    Nonlinear(out, in, w, g::Function, J_g::Function, id=:my_node)
"""
mutable struct Nonlinear <: SoftFactor
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol, Interface}

    g::Function # Vector function that expresses the output vector as a function of the input vector; reduces to scalar for 1-d
    J_g::Function # Jacobi matrix of g, as a function of the input vector; in the 1-d case this reduces to the first derivative of g

    function Nonlinear(out::Variable, in1::Variable, w::Variable, g::Function, J_g::Function; id=ForneyLab.generateId(Nonlinear))
        self = new(id, Vector{Interface}(3), Dict{Symbol,Interface}(), g, J_g)
        ForneyLab.addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

        return self
    end
end

slug(::Type{Nonlinear}) = "g"

# Average energy functional
function averageEnergy(::Type{Nonlinear}, marg_out::ProbabilityDistribution{Univariate}, marg_in1::ProbabilityDistribution{Univariate}, marg_w::ProbabilityDistribution{Univariate}, g::Function, J_g::Function)
    (A, b) = approximate(unsafeMean(marg_in1), g, J_g)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_w) +
    0.5*unsafeMean(marg_w)*( A^2*unsafeCov(marg_in1) + unsafeCov(marg_out) + (A*unsafeMean(marg_in1) + b - unsafeMean(marg_out))^2 )
end

function averageEnergy(::Type{Nonlinear}, marg_out::ProbabilityDistribution{Multivariate}, marg_in1::ProbabilityDistribution{Multivariate}, marg_w::ProbabilityDistribution{MatrixVariate}, g::Function, J_g::Function)
    (A, b) = approximate(unsafeMean(marg_in1), g, J_g)

    0.5*dims(marg_out)*log(2*pi) -
    0.5*unsafeDetLogMean(marg_w) +
    0.5*trace( unsafeMean(marg_w)*( A*unsafeCov(marg_in1)*A' + unsafeCov(marg_out) + (A*unsafeMean(marg_in1) + b - unsafeMean(marg_out))*(A*unsafeMean(marg_in1) + b - unsafeMean(marg_out))' ))
end

function averageEnergy(::Type{Nonlinear}, marg_out_in1::ProbabilityDistribution{Multivariate, Gaussian}, marg_prec::ProbabilityDistribution{Univariate}, g::Function, J_g::Function)
    ensureParameters!(marg_out_in1, (:m, :v))

    V = marg_out_in1.params[:v]
    m = marg_out_in1.params[:m]
    (A, b) = approximate(m[2], g, J_g)

    0.5*log(2*pi) -
    0.5*unsafeLogMean(marg_prec) +
    0.5*unsafeMean(marg_prec)*( A^2*V[2,2] - A*(V[2,1] + V[1,2]) + V[1,1] + (A*m[2] + b - m[1])^2 )
end

function averageEnergy(::Type{Nonlinear}, marg_out_in1::ProbabilityDistribution{Multivariate, Gaussian}, marg_prec::ProbabilityDistribution{MatrixVariate}, g::Function, J_g::Function)
    ensureParameters!(marg_out_in1, (:m, :v))

    V = marg_out_in1.params[:v]
    m = marg_out_in1.params[:m]
    d = Int64(dims(marg_out_in1)/2)
    (A, b) = approximate(m[d+1:end], g, J_g)

    0.5*d*log(2*pi) -
    0.5*unsafeDetLogMean(marg_prec) +
    0.5*trace( unsafeMean(marg_prec)*( A*V[d+1:end,d+1:end]*A' - A*V[d+1:end,1:d] - V[1:d,d+1:end]*A' + V[1:d,1:d] + (A*m[d+1:end] + b - m[1:d])*(A*m[d+1:end] + b - m[1:d])' ) )
end