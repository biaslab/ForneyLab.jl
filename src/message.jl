export Message

type Message{T<:MessagePayload}
    # Message has a payload, which must be a subtype of MessagePayload.
    payload::T
end
Message(payload::MessagePayload) = Message{typeof(payload)}(deepcopy(payload))
Message() = Message(1.0)
==(msg1::Message, msg2::Message) = msg1.payload == msg2.payload
show(io::IO, message::Message) = println(io, "$(typeof(message)) with payload $(message.payload)")