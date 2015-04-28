export  Algorithm, current_algorithm

type Algorithm
    execute::Function
    fields::Dict{Symbol, Any}

    function Algorithm(execute::Function, fields::Dict)
    	self = new(execute, fields)
    	global current_algorithm = self
    	return self
	end
end

global current_algorithm = nothing

function Algorithm(schedule::Schedule, graph::FactorGraph=currentGraph())
    # Constructs an algorithm that executes a predifined schedule
    exec(fields) = execute(fields[:schedule])
    return Algorithm(exec, {:schedule => schedule})
end
