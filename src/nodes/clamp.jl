############################################
# ClampNode
############################################
# Description:
#   ClampNode is used to clamp a parameter to a fixed value.
#   It has only one interface: out.
#   This type is not exported and for internal use only.
############################################

type ClampNode <: Node
    interfaces::Array{Interface, 1}
    out::Interface

    function ClampNode()
        self = new(Array(Interface,1))
        self.interfaces[1] = Interface(self)
        self.out = self.interfaces[1]
        return self
    end
end

function ClampNode(message::Message)
    self = ClampNode()
    self.out.message = message
    return self
end
