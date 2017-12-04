export
ProbabilityDistribution,
Univariate,
Multivariate,
MatrixVariate,
PointMass,
@~,
mean,
var,
cov,
differentialEntropy,
averageEnergy,
==,
vague,
sample,
dims

abstract VariateType
abstract Univariate <: VariateType
abstract Multivariate <: VariateType
abstract MatrixVariate <: VariateType

"""Encodes a probability distribution as a FactorNode of type `family` with fixed interfaces"""
immutable ProbabilityDistribution{var_type<:VariateType, family<:FactorNode}
    params::Dict
end

"""Extract VariateType from dist"""
variateType{V<:VariateType, F<:FactorNode}(dist::ProbabilityDistribution{V, F}) = V

show(io::IO, dist::ProbabilityDistribution) = println(io, format(dist))

matches{Pa<:ProbabilityDistribution, Pb<:ProbabilityDistribution}(Ta::Type{Pa}, Tb::Type{Pb}) = (Pa<:Pb)
matches{T<:ProbabilityDistribution}(::Type{Void}, ::Type{T}) = false

mean(dist::ProbabilityDistribution) = isProper(dist) ? unsafeMean(dist) : error("mean($(dist)) is undefined because the distribution is improper.")
var(dist::ProbabilityDistribution) = isProper(dist) ? unsafeVar(dist) : error("var($(dist)) is undefined because the distribution is improper.")
cov(dist::ProbabilityDistribution) = isProper(dist) ? unsafeCov(dist) : error("cov($(dist)) is undefined because the distribution is improper.")

"""
PointMass is an abstract type used to describe point mass distributions.
It never occurs in a FactorGraph, but it is used as a probability distribution type.
"""
abstract PointMass <: DeltaFactor

slug(::Type{PointMass}) = "δ"

format{V<:VariateType}(dist::ProbabilityDistribution{V, PointMass}) = "$(slug(PointMass))(m=$(format(dist.params[:m])))"

dims(dist::ProbabilityDistribution{Univariate, PointMass}) = 1
dims(dist::ProbabilityDistribution{Multivariate, PointMass}) = length(dist.params[:m])
dims(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = size(dist.params[:m])

# PointMass distribution constructors
ProbabilityDistribution(::Type{Univariate}, ::Type{PointMass}; m::Number=1.0) = ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>m))
ProbabilityDistribution(::Type{PointMass}; m::Number=1.0) = ProbabilityDistribution{Univariate, PointMass}(Dict(:m=>m))
ProbabilityDistribution(::Type{Multivariate}, ::Type{PointMass}; m::Vector=[1.0]) = ProbabilityDistribution{Multivariate, PointMass}(Dict(:m=>m))
ProbabilityDistribution(::Type{MatrixVariate}, ::Type{PointMass}; m::AbstractMatrix=[1.0].') = ProbabilityDistribution{MatrixVariate, PointMass}(Dict(:m=>m))

unsafeMean{T<:VariateType}(dist::ProbabilityDistribution{T, PointMass}) = deepcopy(dist.params[:m])

function unsafeMeanVector(dist::ProbabilityDistribution{Univariate, PointMass}) # For observed bernoulli
    p = dist.params[:m]
    if length(p) == 1
        return [first(p), 1 - first(p)]
    else
        return p
    end
end

unsafeInverseMean(dist::ProbabilityDistribution{Univariate, PointMass}) = 1.0/dist.params[:m]
unsafeInverseMean(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = cholinv(dist.params[:m])

unsafeLogMean(dist::ProbabilityDistribution{Univariate, PointMass}) = log(dist.params[:m])
unsafeLogMean(dist::ProbabilityDistribution{Multivariate, PointMass}) = log(dist.params[:m])
unsafeDetLogMean(dist::ProbabilityDistribution{MatrixVariate, PointMass}) = log(det(dist.params[:m]))

unsafeMirroredLogMean(dist::ProbabilityDistribution{Univariate, PointMass}) = log(1.0 - dist.params[:m])

unsafeVar(dist::ProbabilityDistribution{Univariate, PointMass}) = 0.0
unsafeVar(dist::ProbabilityDistribution{Multivariate, PointMass}) = zeros(dims(dist))

unsafeCov(dist::ProbabilityDistribution{Univariate, PointMass}) = 0.0
unsafeCov(dist::ProbabilityDistribution{Multivariate, PointMass}) = zeros(dims(dist), dims(dist))

isProper{T<:VariateType}(::ProbabilityDistribution{T, PointMass}) = true

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

    # Build the expression for constructing the node id, either an id is passed; or we generate our own
    if isa(dist_expr.args[end], Expr) && (dist_expr.args[end].head == :kw) && (dist_expr.args[end].args[1] == :id)
        node_id_expr = dist_expr.args[end].args[2]
    else
        node_id_expr = parse("ForneyLab.generateId(Variable)")
    end

    expr = parse("""
                begin
                # Use existing object if it exists, otherwise create a new Variable
                $(variable_expr) = try $(variable_expr) catch Variable(id=ForneyLab.pack($(node_id_expr))) end

                # Create new variable if:
                #   - the existing object is not a Variable
                #   - the existing object is a Variable from another FactorGraph
                if (!isa($(variable_expr), Variable)
                    || !haskey(currentGraph().variables, $(variable_expr).id)
                    || !is(currentGraph().variables[$(variable_expr).id], $(variable_expr)))

                    $(variable_expr) = Variable(id=ForneyLab.pack($(node_id_expr)))
                end

                $(dist_expr)
                $(variable_expr)
                end
            """)
    return esc(expr)
end

"""
This function ensures the argument expression is evaluated at runtime, allowing access to local variables
"""
pack(expr) = expr

=={fam_t, fam_u, var_t, var_u}(t::ProbabilityDistribution{var_t, fam_t}, u::ProbabilityDistribution{var_u, fam_u}) = (fam_t == fam_u) && (var_t == var_u) && (t.params == u.params)

function gaussianQuadrature(h::Function; m::Float64=0.0, v::Float64=1.0)
    abscissas = [-7.12581, -6.4095, -5.81223, -5.27555, -4.77716, -4.30555, -3.85376, -3.41717, -2.99249, -2.57725, -2.1695, -1.76765, -1.37038, -0.9765, -0.584979, -0.194841, 0.194841, 0.584979, 0.9765, 1.37038, 1.76765, 2.1695, 2.57725, 2.99249, 3.41717, 3.85376, 4.30555, 4.77716, 5.27555, 5.81223, 6.4095, 7.12581]
    weights = [7.31068e-23, 9.23174e-19, 1.19734e-15, 4.21501e-13, 5.93329e-11, 4.09883e-9, 1.57417e-7, 3.65059e-6, 5.41658e-5, 0.000536268, 0.00365489, 0.0175534, 0.0604581, 0.15127, 0.277458, 0.375238, 0.375238, 0.277458, 0.15127, 0.0604581, 0.0175534, 0.00365489, 0.000536268, 5.41658e-5, 3.65059e-6, 1.57417e-7, 4.09883e-9, 5.93329e-11, 4.21501e-13, 1.19734e-15, 9.23174e-19, 7.31068e-23]

    1/sqrt(pi) * sum([weights[i]*h(sqrt(2*v)*abscissas[i] + m) for i=1:32])
end

"""isValid: return true if the parameter field exists and (the first element of) the parameter is not NaN"""
isValid(dist::ProbabilityDistribution, field::Symbol) = ( haskey(dist.params, field) && !isnan(dist.params[field][1]) )

function invalidate!(dist::ProbabilityDistribution, field::Symbol)
    if haskey(dist.params, field)
        dist.params[field][1] = NaN
    end
    return dist
end