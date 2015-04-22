export  Algorithm

type Algorithm
    execute::Function
    fields::Dict{Symbol, Any}
end

global current_algorithm = nothing

function Algorithm(schedule::Schedule, graph::FactorGraph=currentGraph())
    # Constructs an algorithm that executes a predifined schedule
    exec(fields) = execute(fields[:schedule])
    return Algorithm(exec, {:schedule => schedule})
end
