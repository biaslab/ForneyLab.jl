export
    GaussianMessage,
    GeneralMessage

############################################
# GaussianMessage
############################################
# Description:
#   Encodes a Gaussian PDF with covariance V
#   and mean m.
############################################
type GaussianMessage <: Message
    V::Array{Float64}
    m::Array{Float64}
end
GaussianMessage() = GaussianMessage([1.0], [0.0])

############################################
# GeneralMessage
############################################
# Description:
#   Simply holds an arbitrary object.
#   Useful for example for passing parameters.
############################################
type GeneralMessage{T} <: Message
    value::T
end
GeneralMessage() = GeneralMessage(1.0)