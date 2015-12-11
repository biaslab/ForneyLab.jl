export AbstractMessage, EmptyMessage, Message

import Base.isempty

abstract AbstractMessage

type EmptyMessage <: AbstractMessage end

type Message{T<:ProbabilityDistribution} <: AbstractMessage
    # Message has a payload, which must be a ProbabilityDistribution.
    payload::T
end

==(msg1::Message, msg2::Message) = (msg1.payload == msg2.payload)
==(msg1::EmptyMessage, msg2::Message) = false
==(msg1::Message, msg2::EmptyMessage) = false
==(msg1::EmptyMessage, msg2::EmptyMessage) = true

show(io::IO, message::EmptyMessage) = println(io, "Empty Message")
show(io::IO, message::Message) = println(io, "$(typeof(message)) with payload $(message.payload)")
