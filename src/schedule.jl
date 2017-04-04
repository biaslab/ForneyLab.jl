"""
A Schedule is a list of `(interface, msg_update_rule)` tuples, and defines
a message passing schedule.
The `msg_update_rule <: MessageUpdateRule` defines the rule that is used
to calculate the message coming out of `interface`.
"""
typealias Schedule Vector{Tuple{Interface, DataType}}