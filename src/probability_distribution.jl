export
ProbabilityDistribution,
PointMass,
@~,
mean,
var,
cov,
differentialEntropy,
averageEnergy,
==,
vague

# TODO: correctness of distribution parameters is not enforced
"""Encodes a probability distribution as a FactorNode of type `family` with fixed interfaces"""
immutable ProbabilityDistribution{family<:FactorNode}
    params::Dict
end

ProbabilityDistribution{T<:SoftFactor}(family::Type{T}; kwargs...) = ProbabilityDistribution{family}(Dict(kwargs))

mean(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
var(dist::ProbabilityDistribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::ProbabilityDistribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")

"""
PointMass is an abstract type used to describe point mass distributions.
It never occurs in a FactorGraph, but it is used as a probability distribution type.
"""
abstract PointMass <: DeltaFactor

ProbabilityDistribution(family::Type{PointMass}; kwargs...) = ProbabilityDistribution{family}(Dict(kwargs))

vague(::Type{ProbabilityDistribution{PointMass}}) = ProbabilityDistribution(PointMass, m=1.0)

unsafeMean(dist::ProbabilityDistribution{PointMass}) = dist.params[:m]

unsafeInverseMean(dist::ProbabilityDistribution{PointMass}) = 1.0/dist.params[:m]

unsafeLogMean(dist::ProbabilityDistribution{PointMass}) = log(dist.params[:m])

unsafeMirroredLogMean(dist::ProbabilityDistribution{PointMass}) = log(1.0 - dist.params[:m])

unsafeVar(dist::ProbabilityDistribution{PointMass}) = 0.0

unsafeCov(dist::ProbabilityDistribution{PointMass}) = 0.0

mean(dist::ProbabilityDistribution{PointMass}) = unsafeMean(dist)

isProper(::ProbabilityDistribution{PointMass}) = true

"""
`x ~ GaussianMeanVariance(constant(0.0), constant(1.0), id=:some_optional_id)` is a shorthand notation
for `x = Variable(); GaussianMeanVariance(x, constant(0.0), constant(1.0))`
"""
macro ~(variable_expr::Any, dist_expr::Expr)
    # Sanity checks
    (dist_expr.head == :call) || error("Incorrect use of ~ operator.")
    (eval(dist_expr.args[1]) <: SoftFactor) || error("~ operator should be followed by subtype of SoftFactor.")

    # Build FactorNode constructor call
    if isa(dist_expr.args[2], Expr) && (dist_expr.args[2].head == :parameters)
        dist_expr.args = vcat(dist_expr.args[1:2], [variable_expr], dist_expr.args[3:end])
    else
        dist_expr.args = vcat([dist_expr.args[1]; variable_expr], dist_expr.args[2:end])
    end

    expr = parse("""
                begin
                # Use existing object if it exists, otherwise create a new Variable
                $(variable_expr) = try $(variable_expr) catch Variable() end

                # Create new variable if:
                #   - the existing object is not a Variable
                #   - the existing object is a Variable from another FactorGraph
                if (!isa($(variable_expr), Variable)
                    || !haskey(currentGraph().variables, $(variable_expr).id)
                    || !is(currentGraph().variables[$(variable_expr).id], $(variable_expr)))

                    $(variable_expr) = Variable()
                end

                $(dist_expr)
                $(variable_expr)
                end
            """)
    return esc(expr)
end

=={T, U}(t::ProbabilityDistribution{T}, u::ProbabilityDistribution{U}) = (T == U) && (t.params == u.params)

function gaussianQuadrature(h::Function; m::Float64=0.0, v::Float64=1.0)
    abscissas = [-7.12581, -6.4095, -5.81223, -5.27555, -4.77716, -4.30555, -3.85376, -3.41717, -2.99249, -2.57725, -2.1695, -1.76765, -1.37038, -0.9765, -0.584979, -0.194841, 0.194841, 0.584979, 0.9765, 1.37038, 1.76765, 2.1695, 2.57725, 2.99249, 3.41717, 3.85376, 4.30555, 4.77716, 5.27555, 5.81223, 6.4095, 7.12581]
    weights = [7.31068e-23, 9.23174e-19, 1.19734e-15, 4.21501e-13, 5.93329e-11, 4.09883e-9, 1.57417e-7, 3.65059e-6, 5.41658e-5, 0.000536268, 0.00365489, 0.0175534, 0.0604581, 0.15127, 0.277458, 0.375238, 0.375238, 0.277458, 0.15127, 0.0604581, 0.0175534, 0.00365489, 0.000536268, 5.41658e-5, 3.65059e-6, 1.57417e-7, 4.09883e-9, 5.93329e-11, 4.21501e-13, 1.19734e-15, 9.23174e-19, 7.31068e-23]

    1/sqrt(pi) * sum([weights[i]*h(sqrt(2*v)*abscissas[i] + m) for i=1:32])
end