@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{LogNormal},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPExponentialOutNG)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{LogNormal}, Nothing),
                :name          => SPExponentialIn1LN)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPExponentialOutNP)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPExponentialIn1PN)