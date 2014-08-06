############################################
# ClampNode
############################################
# Description:
#   ClampNode is used to clamp a parameter to a fixed value.
#   It has only one out interface.
#   This type is not exported and for internal use only.
############################################

type ClampNode <: Node
    out::Interface

    function ClampNode()
        self = new()
        self.out = Interface(self)
        return self
    end
end

function ClampNode(message::Message)
    self = ClampNode()
    self.out.message = message
    return self
end
