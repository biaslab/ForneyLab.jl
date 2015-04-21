    # collectInbounds(...) as sum-product specific function.
    # for j = 1:length(node.interfaces)
    #     interface = node.interfaces[j]
    #     if interface == outbound_interface
    #         outbound_interface_id = j
    #         push!(inbound_array, nothing) # This interface is outbound, push "nothing"
    #     else
    #         # Inbound message or marginal is required, push the required message/marginal to inbound_array
    #         pushRequiredInbound!(algorithm, inbound_array, node, interface, outbound_interface)
    #     end
    # end
