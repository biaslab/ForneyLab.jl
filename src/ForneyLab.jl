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
ScalarParameterMessage() = ScalarParameterMessage(1.0)

#
# Nodes
#
abstract Node

# An interface is a half-edge that connects to another interface to form a complete edge.
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

# Combine two interfaces to form a complete edge.
type Edge
    tail::Union(Interface, Nothing)
    head::Union(Interface, Nothing)

    function Edge(tail::Interface, head::Interface)
        if  typeof(tail) == Nothing ||
            typeof(head) == Nothing ||
            typeof(head.message) == typeof(tail.message)

            tail.partner = head
            head.partner = tail
            new(tail, head)
        else
            error("Head and tail message types do not match")
        end
    end
end

#
# calculatemessage calculates the outgoing message on the sink interface (half-edge).
# Every combination of node and message types will get an individual signature for the function call.
#
function calculatemessage(  sink::Interface,
                            node::MultiplicationNode,
                            multiplierMessage::ScalarParameterMessage,
                            sourceMessage::GaussianMessage)
    # Multiplication of a Gaussian message
    sink_V = multiplierMessage.value' * sourceMessage.V * multiplierMessage.value
    sink_m = multiplierMessage.value * sourceMessage.m
    sink.message = GaussianMessage(sink_V, sink_m)
end

# call for calculatemessage that calculates the incoming messages by recursion.
function calculatemessage(  sink::Interface,
                            node::MultiplicationNode)
    # TODO: Collect all incoming messages
    calculatemessage(sink, node, ScalarParameterMessage(4.0), GaussianMessage())
end

calculatemessage(sink::Interface) = calculatemessage(sink, sink.node)

calculateforwardmessage(edge::Edge) = calculatemessage(edge.tail)
calculatebackwardmessage(edge::Edge) = calculatemessage(edge.head)

end # module ForneyLab