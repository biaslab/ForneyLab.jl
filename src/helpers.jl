export getArgumentValue, isApproxEqual

# ensureMatrix: ensure that the input is a 2D array or nothing
ensureMatrix{T<:Number}(arr::Array{T, 2}) = arr
ensureMatrix{T<:Number}(arr::Array{T, 1}) = reshape(arr, 1, 1)
ensureMatrix(n::Nothing) = nothing

# getArgumentValue: get value of specific keyword argument
function getArgumentValue(args::Array, key::Symbol)
    for (k, v) in args
        if k==key
            return v
        end
    end
    return false
end

# isApproxEqual: check approximate equality
isApproxEqual(arg1, arg2) = maximum(abs(arg1-arg2)) < 1.0e-12

# isRoundedPosDef: is input matrix positive definite? Round to prevent fp precision problems that isposdef() suffers from.
isRoundedPosDef{T<:FloatingPoint}(arr::Array{T, 2}) = ishermitian(round(arr, 12)) && isposdef(arr, 'L')

