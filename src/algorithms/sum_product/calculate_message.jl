function calculateMessage!(outbound_interface::Interface)
    # Calculate the outbound message on a specific interface by generating a schedule and executing it.
    # The resulting message is stored in the specified interface and returned.

    # Lock graph structure
    scheme = InferenceScheme()

    # Generate a message passing schedule
    printVerbose("Auto-generating message passing schedule...\n")
    schedule = generateSchedule!(outbound_interface, scheme)
    if verbose show(schedule) end

    # Execute the schedule
    printVerbose("\nExecuting above schedule...\n")
    execute(schedule, scheme)
    printVerbose("\ncalculateMessage!() done.")

    return outbound_interface.message
end

# Calculate forward/backward messages on an Edge
calculateForwardMessage!(edge::Edge) = calculateMessage!(edge.tail)
calculateBackwardMessage!(edge::Edge) = calculateMessage!(edge.head)