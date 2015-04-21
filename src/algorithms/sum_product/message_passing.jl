function collectInbounds(outbound_interface::Interface)
    # Sum-product specific method to collect all required inbound messages in an array.
    # This array is used to call the node update function (sumProduct!)
    # outbound_interface: the interface on which the outbound message will be updated
    # Returns: (outbound interface id, array of inbound messages)
    
    outbound_interface_id = 0
    inbounds = Array(Any, 0)
    for j = 1:length(outbound_interface.node.interfaces)
        interface = node.interfaces[j]
        if is(interface, outbound_interface)
            # We don't need the inbound message on the outbound interface
            outbound_interface_id = j
            push!(inbounds, nothing) # This interface is outbound, push "nothing"
        else
            try 
                push!(inbounds, interface.partner.message) 
            catch 
                error("Cannot collect inbound message on $(interface). Make sure there is an inbound message present at this interface.") 
            end
        end
    end

    return (outbound_interface_id, inbounds)
end