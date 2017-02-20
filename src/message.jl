export Message

"""
A MessageCalculationRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageCalculationRule)` to list the available rules.
"""
abstract MessageCalculationRule

abstract SumProductRule <: MessageCalculationRule
abstract VariationalRule <: MessageCalculationRule
abstract ExpectationRule <: MessageCalculationRule

"""
TODO: description
"""
type Message{rule<:MessageCalculationRule}
    interface::Interface
    rule_id::Symbol
    outbound_family::DataType

    Message{rule<:MessageCalculationRule}(interface::Interface) = new{rule}(interface)
end