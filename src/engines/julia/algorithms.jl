export sumProductAlgorithm, variationalAlgorithm, expectationPropagationAlgorithm, variationalExpectationPropagationAlgorithm

"""
Create a sum-product algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function sumProductAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    schedule = sumProductSchedule(variables)
    algo = messagePassingAlgorithm(schedule, variables, file=file, name=name)

    return algo
end
sumProductAlgorithm(variable::Variable; file::String="", name::String="") = sumProductAlgorithm([variable], file=file, name=name)

"""
Create a variational algorithm to infer marginals over a recognition distribution, and compile it to Julia code
"""
function variationalAlgorithm(q_factor::RecognitionFactor; file::String="", name::String="")
    q_schedule = variationalSchedule(q_factor)
    algo = messagePassingAlgorithm(q_schedule, collect(q_factor.variables), file=file, name=name)

    return algo
end

"""
Create a variational algorithm to infer marginals over the recognition distribution that includes `variables`, and compile it to Julia code
"""
function variationalAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    q_factor = RecognitionFactor(variables)

    return variationalAlgorithm(q_factor, file=file, name=name)
end
variationalAlgorithm(variable::Variable; file::String="", name::String="") = variationalAlgorithm([variable], file=file, name=name)

"""
Create an expectation propagation algorithm to infer marginals over `variables`, and compile it to Julia code
"""
function expectationPropagationAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    schedule = expectationPropagationSchedule(variables)
    algo = messagePassingAlgorithm(schedule, variables, file=file, name=name)

    return algo
end
expectationPropagationAlgorithm(variable::Variable; file::String="", name::String="") = expectationPropagationAlgorithm([variable], file=file, name=name)

"""
Create a variational EP algorithm to infer marginals over the recognition distribution that includes `variables`, and compile it to Julia code
"""
function variationalExpectationPropagationAlgorithm(variables::Vector{Variable}; file::String="", name::String="")
    q_factor = RecognitionFactor(variables)
    q_schedule = variationalExpectationPropagationSchedule(q_factor)
    algo = messagePassingAlgorithm(q_schedule, collect(q_factor.variables), file=file, name=name)

    return algo
end
variationalExpectationPropagationAlgorithm(variable::Variable; file::String="", name::String="") = variationalExpectationPropagationAlgorithm([variable], file=file, name=name)