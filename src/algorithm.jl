export  Algorithm

type Algorithm
    execute::Function
    fields::Dict{Symbol, Any}

    function Algorithm(execute::Function, fields::Dict = Dict{Symbol,Any}())
        self = new(execute, fields)
        global current_algorithm = self
        return self
    end
end

function Algorithm(schedule::Schedule, graph::FactorGraph=currentGraph())
    # Constructs an algorithm that executes a predifined schedule
    exec(fields) = execute(fields[:schedule])
    return Algorithm(exec, Dict{Any,Any}(:schedule => schedule))
end

function show(io::IO, algo::Algorithm)
    println(io, "Algorithm with fields:")
    for (key, val) in algo.fields
        println(io, " $(key)::$(typeof(val))")
    end
    println(io, "\nUse algorithm.fields[:field] to inspect field values.")
end

Base.deepcopy(::Algorithm) = error("deepcopy(::Algorithm) is not possible. You should construct a new Algorithm.")
