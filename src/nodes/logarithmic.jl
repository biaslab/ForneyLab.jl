############################################
# LogarithmicNode
############################################
# Description:
#   Maps a Gaussian to a gamma distribution.
#
#    in1        out
#   ----->[log]----->
#
#   Example:
#       LogarithmicNode(; name="my_node")
#
# Interface ids, (names) and supported message types:
#   1. (in1):
#       Message{GaussianDistribution}
#   2. (out):
#       Message{GammaDistribution}
############################################