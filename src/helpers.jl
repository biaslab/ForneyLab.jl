export isApproxEqual, huge, tiny, KLpq

# ensureMatrix: ensure that the input is a 2D array or nothing
ensureMatrix{T<:Number}(arr::Array{T, 2}) = arr
ensureMatrix{T<:Number}(arr::Array{T, 1}) = reshape(arr, 1, 1)
ensureMatrix(n::Nothing) = nothing

# Define what is considered big and small, mostly used for setting vague messages.
huge(::Type{Float64}) = 1e12
huge() = huge(Float64)
tiny(::Type{Float64}) = 1./huge(Float64)
tiny() = tiny(Float64)

# isApproxEqual: check approximate equality
isApproxEqual(arg1, arg2) = maximum(abs(arg1-arg2)) < tiny()

# isRoundedPosDef: is input matrix positive definite? Round to prevent fp precision problems that isposdef() suffers from.
isRoundedPosDef{T<:FloatingPoint}(arr::Array{T, 2}) = ishermitian(round(arr, int(log(10, huge())))) && isposdef(arr, :L)

function viewFile(filename::String)
    # Open a file with the application associated with the file type
    @windows? run(`cmd /c start $filename`) : (@osx? run(`open $filename`) : (@linux? run(`xdg-open $filename`) : error("Cannot find an application for $filename")))
end

global unnamed_counter = 0
function unnamedStr()
    global unnamed_counter::Int64 += 1
    return "unnamed$(unnamed_counter)"
end

function KLpq(x::Vector{Float64}, p::Vector{Float64}, q::Vector{Float64})
    # Calculate the KL divergence of p and q
    length(x) == length(p) == length(q) || error("Lengths of x, p and q must match")
    dx = diff(x)
    return -sum(p[1:end-1].*log(q[1:end-1]./p[1:end-1]).*dx)
end
