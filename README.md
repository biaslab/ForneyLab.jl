ForneyLab.jl
============

Forney-style Factor Graph toolbox in Julia.
**This software is still very experimental**.

Installation
============
Add the ForneyLab package to your Julia installation:
```jl
Pkg.clone("https://github.com/spsbrats/ForneyLab.jl.git")
```
To update an already installed package, use:
```jl
Pkg.update()
```
Usage
=====
Import ForneyLab:
```jl
using ForneyLab
```
Once imported, one can create nodes and edges to build a factor graph. There are multiple types of nodes. One can use built-in node types, or one can define custom nodes. A node has one or more interfaces, which can be used to send/receive messages. An edge connects two interfaces of different nodes. Every interface can handle one or more message types. ForneyLab comes with a set of built-in message types, but you can also define your own.

The demo directory contains some illustrative demos to get you started.

Example
-------
TODO: add description.
```jl
using ForneyLab
# Build a graph
c1 = ConstantNode(GeneralMessage(3.0))
c2 = ConstantNode(GeneralMessage(2.0))
x = MultiplicationNode(name="MyMultiplier")
e1 = Edge(c1.out, x.in1)
e2 = Edge(c2.out, x.in2)
# Calculate the outbound message on x.out, the output of the multiplication node.
# This call will recursively calculate all required inbound messages.
calculateMessage!(x.out)
# Print the calculated message, which is stored in the interface.
print(x.out.message)
# We can also request directional messages on edges instead of an interface.
# I.e. for calculating the message c1.interface -> x.in1, we can use:
calculateForwardMessage!(e1)
# The calculated message is again stored in the sending interface:
print(e1.tail.message)
```
Message types
=============
A message type is always a subtype of `Message`, and its name ends in "Message". Built-in message types are defined in "messages.jl". It is easy to add custom message types. For example, one can add a message to hold a binomial distribution:
```jl
type BinomialMessage <: Message
    p::Float64
end
BinomialMessage() = BinomialMessage(0.5)
```
To create a message of type BinomialMessage, use `m = BinomialMessage(0.2)` or just `m = BinomialMessage()` to init p at 0.5. There should always be a message constructor that doesn't require any arguments, for testing purposes.

Node types
==========
A node type is always a subtype of `Node`, and its name ends in "Node". Built-in node types are defined in "nodes/*.jl". A node type at least has fields "interfaces" and "name". For example, an addition node with 3 interfaces can be defined by:
```jl
type AdditionNode <: Node
    interfaces::Array{Interface,1}
    name::ASCIIString
    in1::Interface
    in2::Interface
    out::Interface
    function AdditionNode(;args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end
        self = new(Array(Interface, 3), name)
        # Create interfaces
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Init named interface handles
        self.in1    = self.interfaces[1]
        self.in2    = self.interfaces[2]
        self.out    = self.interfaces[3]
        return self
    end
end
```
Every interface has a unique id, given by its index in the `interfaces` array.
Apart from the node type definition, one also has to define one or more methods for calculating the outbound messages. For calculating messages, function `updateNodeMessage!()` is used. Multiple methods of this function can be defined if one wants separate implementations for different message types. The `updateNodeMessage!()` method for `AdditionNode` could be defined as:
```jl
function updateNodeMessage!{T<:Union(GaussianMessage, GeneralMessage)}
                            (outbound_interface_id::Int,
                             node::AdditionNode,
                             inbound_messages::Array{T, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end

    # TODO: CALCULATE THE OUTBOUND MESSAGE HERE,
    # based on the outbound_interface_id and the inbound_messages.
    # The calculated message is saved in node.interfaces[outbound_interface_id].message and should also be returned.

    return node.interfaces[outbound_interface_id].message = GaussianMessage()
end
```
This method will be used when calculating outbound messages of a AdditionNode, when the inbound messages are of types GaussianMessage and/or GeneralMessage. For different message types, extra methods can be defined.
