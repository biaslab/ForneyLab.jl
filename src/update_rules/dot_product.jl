#             β
#             | ↓ N{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        →        →
#       δ{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{Gaussian}),
                :name          => SPDotProductOutNPG)


#             β
#             | ↓ δ{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        →        →
#       N{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPDotProductOutNGP)


#             β
#             | ↑ N{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        →        ←
#       δ{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Nothing),
                :name          => SPDotProductIn2GPN)


#             β
#             | ↓ δ{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        ←        ←
#       N{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPDotProductIn1GNP)