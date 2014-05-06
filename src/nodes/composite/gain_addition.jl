############################################
# GainAdditionCompositeNode
############################################
# Description:
#   Composite gain-addition node.
#
#            | in1
#        ____|____
#        |   |   |
#        |  [A]  |
#        |   |   |
#    in2 |   v   | out
#   -----|->[+]--|---->
#        |_______|
#
#
#   out = A*in1 + in2
#   Example:
#       GainAdditionCompositeNode([1.0]; name="my_node")
#   Gain A is fixed and defined through the constructor.
#
# Interface ids, (names) and supported message types:
#   1. in1:
#       GaussianMessage
#   2. in2:
#       GaussianMessage
#   3. out:
#       GaussianMessage
#
############################################
