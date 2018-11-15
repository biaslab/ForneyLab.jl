import Base: -
export -

function -(in1::Variable, in2::Variable)
    out = Variable()
    Addition(in1, out, in2)
    return out
end

function -(in1::Variable, in2)
    @ensureVariables(in2)
    in1 - in2
end

-(in1, in2::Variable) = -(in2, in1)
