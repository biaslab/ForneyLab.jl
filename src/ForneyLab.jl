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

#
# Nodes
#
abstract Node

type Interface
    node::Node
    partner::Interface
end

# Mulitiply with parameter a: dst = a * src
type MultiplicationNode <: Node
    a::Interface
    src::Interface
    dst::Interface
    function MultiplicationNode()
        self = new()
        a = Interface(self)
        src = Interface(self)
        dst = Interface(self)
    end
end

end # module ForneyLab