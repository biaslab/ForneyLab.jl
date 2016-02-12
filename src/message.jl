export AbstractMessage, Message

import Base.isempty

"""
A Message sits on an interface and has a payload, which must be a ProbabilityDistribution.
"""
type Message{T<:ProbabilityDistribution}
    payload::T
end

==(msg1::Message, msg2::Message) = (msg1.payload == msg2.payload)

show(io::IO, message::Message) = println(io, "$(typeof(message)) with payload $(message.payload)")
