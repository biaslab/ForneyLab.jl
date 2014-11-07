export Message

type Message{T<:ProbabilityDistribution}
    # Message has a payload, which must be a subtype of ProbabilityDistribution.
    payload::T
end

==(msg1::Message, msg2::Message) = msg1.payload == msg2.payload
show(io::IO, message::Message) = println(io, "$(typeof(message)) with payload $(message.payload)")