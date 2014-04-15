export GaussianMessage, ScalarParameterMessage

type GaussianMessage <: Message
    # Encodes a Gaussian PDF with covariance V and mean m
    V::Array{Float64}
    m::Array{Float64}
end
GaussianMessage() = GaussianMessage([1.0], [0.0])

type ScalarParameterMessage{T} <: Message
    # Simply holds a variable
    value::T
end
ScalarParameterMessage() = ScalarParameterMessage(1.0)