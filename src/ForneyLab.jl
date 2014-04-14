module ForneyLab

#
# Messages
#
abstract Message

type GaussianMessage <: Message
    V::Array{Float64}
    m::Array{Float64}
end
GaussianMessage() = GaussianMessage([1.0], [0.0])

type ScalarParameterMessage{T} <: Message
    value::T
end

#
# Nodes
#
abstract Node

type Interface
    node::Node
    partner::Union(Interface, Nothing)
    message::Union(Message, Nothing)
end
Interface(node::Node) = Interface(node, nothing, nothing)

# Mulitiply source with parameter multiplier: sink = multiplier * source
type MultiplicationNode <: Node
    multiplier::Interface
    source::Interface
    sink::Interface
    function MultiplicationNode()
        self = new()
        self.multiplier = Interface(self)
        self.source = Interface(self)
        self.sink = Interface(self)
        return self
    end
end

function calculatemessage(  sink::Interface,
                            node::MultiplicationNode,
                            multiplierMessage::ScalarParameterMessage,
                            sourceMessage::GaussianMessage)
    sink.message = GaussianMessage()
end

function calculatemessage(  sink::Interface,
                            node::MultiplicationNode)
    # TODO: Collect all incoming messages
    calculatemessage(sink, node, ScalarParameterMessage(4.0), GaussianMessage())
end

calculatemessage(sink::Interface) = calculatemessage(sink, sink.node)

end # module ForneyLab