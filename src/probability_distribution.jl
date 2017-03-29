export
ProbabilityDistribution,
PointMass,
@~,
mean,
var,
cov

"""Encodes a probability distribution as a FactorNode of type `family` with fixed interfaces"""
immutable ProbabilityDistribution{family<:FactorNode}
    parameters::Tuple
end

ProbabilityDistribution{T<:SoftFactor}(family::Type{T}, args...) = ProbabilityDistribution{family}(args)

mean(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
var(dist::ProbabilityDistribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::ProbabilityDistribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")

"""
PointMass is an abstract type used to describe point mass distributions.
It never occurs in a FactorGraph, but it is used as a probability distribution type.
"""
abstract PointMass <: DeltaFactor

ProbabilityDistribution(family::Type{PointMass}, args...) = ProbabilityDistribution{family}(args)

mean(dist::ProbabilityDistribution{PointMass}) = dist.parameters[1]

isProper(::ProbabilityDistribution{PointMass}) = true

"""
`x ~ GaussianMeanVariance(0.0, 1.0, id=:some_optional_id)` is a shorthand notation
for `x = Variable(); GaussianMeanVariance(x, constant(0.0), constant(1.0))`
"""
macro ~(variable_id::Symbol, dist_expr::Expr)
    # Sanity checks
    (dist_expr.head == :call) || error("Incorrect use of ~ operator.")
    (eval(dist_expr.args[1]) <: FactorNode) || error("~ operator should be followed by subtype of FactorNode.")

    # Build FactorNode constructor call
    dist_expr.args = vcat([dist_expr.args[1]; variable_id], dist_expr.args[2:end])

    ex = parse("""
                begin
                $(variable_id) = try $(variable_id) catch Variable(id=:($(variable_id))) end
                if !haskey(currentGraph().variables, :($(variable_id))) || !is(currentGraph().variables[:($(variable_id))], $(variable_id))
                    $(variable_id) = Variable(id=:($(variable_id)))
                end
                $(dist_expr)
                end
            """)
    esc(ex)
end
