ForneyLab.jl
============

Forney-style Factor Graph toolbox in Julia.
*This software is still very experimental*.

Installation
============
Add the ForneyLab package to your Julia installation:

    Pkg.clone("https://github.com/spsbrats/ForneyLab.jl.git")

To update an already installed package, use:

    Pkg.update()

Usage
=====
Import ForneyLab:

    using ForneyLab
    
Once imported, one can create nodes and edges to build a factor graph. There are multiple types of nodes. One can use built-in node types, or one can define custom nodes. A node has one or more interfaces, which can be used to send/receive messages. An edge connects two interfaces of different nodes. Every interface can handle one or more message types. ForneyLab comes with a set of built-in message types, but you can also define your own. 

Example
-------
TODO: add description.

    using ForneyLab
    # Build a graph
    c1 = ConstantNode(GeneralMessage(3.0))
    c2 = ConstantNode(GeneralMessage(2.0))
    x = MultiplicationNode("MyMultiplier")
    e1 = Edge(c1.interface, x.source)
    e2 = Edge(c2.interface, x.multiplier)
    # Calculate the outbound message on x.sink, the output of the multiplication node.
    # This call will recursively calculate all required inbound messages.
    calculatemessage!(x.sink)
    # Print the calculated message, which is stored in the interface.
    print(x.sink.message)
    # We can also request directional messages on edges instead of an interface.
    # I.e. for calculating the message c1.interface -> x.source, we can use:
    calculateforwardmessage!(e1)
    # The calculated message is again stored in the sending interface:
    print(e1.tail.message)
    
Message types
=============
A message type is always a subtype of `Message`, and its name ends in "Message". Built-in message types are defined in "messages.jl". It is easy to add custom message types. For example, one can add a message to hold a binomial distribution:

    type BinomialMessage <: Message
        p::Float64
    end
    BinomialMessage() = BinomialMessage(0.5)
    
To create a message of type BinomialMessage, use `m = BinomialMessage(0.2)` or just `m = BinomialMessage()` to init p at 0.5.

Node types
==========
A node type is always a subtype of `Node`, and its name ends in "Node". Built-in node types are defined in "nodes/*.jl". A node type at least has fields "interfaces" and "name". For example, an addition node with 3 interfaces can be defined by:

    type AdditionNode <: Node
        interfaces::Array{Interface,1}
        name::ASCIIString
        in1::Interface
        in2::Interface
        out::Interface
        function AdditionNode(name::ASCIIString="#undef")
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
    
Every interface has a unique id, given by its index in the `interfaces` array.
Apart from the node type definition, one also has to define one or more methods for calculating the outbound messages. For calculating messages, function `calculatemessage!()` is used. Multiple methods of this function can be defined if one wants separate implementations for different message types. The `calculatemessage!()` function for `AdditionNode` could be defined as:

    function calculatemessage!{T<:Union(GaussianMessage,GeneralMessage)}(
                                interfaceId::Int,
                                node::AdditionNode,
                                inboundMessages::Array{T,1},
                                messageType::DataType)
        # Calculate the output message here, 
        # based on the output interface (interfaceId) and the inboundMessages.
        # messageType is the desired type of the output message.
        # The calculated message is saved in node.interfaces[interfaceId].message
        node.interfaces[interfaceId].message = GaussianMessage()
    end
    
This method will be used when calculating outbound messages of a AdditionNode, when the inbound messages are of types GaussianMessage and/or GeneralMessage. For different message types, other methods can be defined.
