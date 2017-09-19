export
ProbabilityDistribution,
PointMass,
@~,
mean,
var,
cov,
==

# TODO: correctness of distribution parameters is not enforced
# TODO: current use of ProbabilityDistribution{family} leads to very long calling signatures
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

unsafeMean(dist::ProbabilityDistribution{PointMass}) = dist.params[:m]

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