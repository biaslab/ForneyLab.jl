"""
A MessageCalculationRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageCalculationRule)` to list the available rules.
"""
abstract MessageCalculationRule

abstract SumProductRule <: MessageCalculationRule
abstract VariationalRule <: MessageCalculationRule
abstract ExpectationRule <: MessageCalculationRule