export
MessageUpdateRule,
SumProductRule

"""
A MessageCalculationRule specifies how a Message is calculated from the node function and the incoming messages.
Use `subtypes(MessageCalculationRule)` to list the available rules.
"""
abstract MessageUpdateRule

abstract SumProductRule{factor_type} <: MessageUpdateRule