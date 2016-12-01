export Message, MessageCalculationRule
export SumProductRule, VariationalRule, ExpectationRule # <: MessageCalculationRule
export implementation

"""
A Message sits on an interface and has a payload, which must be a ProbabilityDistribution.
"""
type Message{T<:ProbabilityDistribution}
    payload::T
end

==(msg1::Message, msg2::Message) = (msg1.payload == msg2.payload)

show(io::IO, message::Message) = println(io, "$(typeof(message)) with payload $(message.payload)")


"""
A MessageCalculationRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageCalculationRule)` to list the available rules.
"""
abstract MessageCalculationRule

abstract SumProductRule <: MessageCalculationRule
implementation(::Type{SumProductRule}) = sumProductRule!

abstract VariationalRule <: MessageCalculationRule
implementation(::Type{VariationalRule}) = variationalRule!

abstract ExpectationRule <: MessageCalculationRule
implementation(::Type{ExpectationRule}) = expectationRule!