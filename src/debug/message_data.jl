export MessageData, MessageSnapshot

struct MessageSnapshot
    variate::String
    family::String
    params::Dict
end

struct MessageData
    edgeID::Union{String, Symbol}
    type::String
    message::MessageSnapshot
end

MessageSnapshot(message::Message{F, V}) where { F, V }    = MessageSnapshot(string(V), string(F), message.dist.params)

function MessageSnapshot(message::Message{F, V}) where { F <: Function, V } 
    # error("Cannot dump function message [WIP]")
    return MessageSnapshot(string(V), string(F), Dict())
end
