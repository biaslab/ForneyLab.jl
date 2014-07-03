############################################
# LinearCompositeNode
############################################
# Description:
#   Linear variational node for parameter estimation on a linear function.
#
#            a  b      s_N
#            |  |       |
#       |-----------------|
#       |    |  |       | |
#       |    |  -->[N]<-- |
#       |    |      |     |
#       |    v      v     |
#  in1--|-->[a]--->[+]----|-->out
#       |                 |
#       |-----------------|
#
#
#
# Interface ids, (names) and supported message types:
#   1. in1:
#       GaussianMessage
#   2. a:
#       GaussianMessage
#   3. b:
#       GaussianMessage
#   4. s_N:
#       InvertedGammaMessage
#   3. out:
#       GaussianMessage
#
############################################
