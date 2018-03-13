#             β
#             | ↓ N{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        →        →
#       δ{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPDotProductOutVPG)


#             β
#             | ↓ δ{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        →        →
#       N{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPDotProductOutVGP)


#             β
#             | ↑ N{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        →        ←
#       δ{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Void),
                :name          => SPDotProductIn2GPV)


#             β
#             | ↓ δ{d}
#             |
#             V
#    x ----->[⋅]-----> y
#        ←        ←
#       N{d}      N

@sumProductRule(:node_type     => DotProduct,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPDotProductIn1GVP)