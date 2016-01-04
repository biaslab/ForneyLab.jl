export AbstractMessage, Message

import Base.isempty

type Message{T<:ProbabilityDistribution}
    # Message has a payload, which must be a ProbabilityDistribution.
    payload::T
end

==(msg1::Message, msg2::Message) = (msg1.payload == msg2.payload)

show(io::IO, message::Message) = println(io, "$(typeof(message)) with payload $(message.payload)")
