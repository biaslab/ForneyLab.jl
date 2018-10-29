@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{LogNormal},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPExponentialOutVG)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{LogNormal}, Nothing),
                :name          => SPExponentialIn1LV)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPExponentialOutVP)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPExponentialIn1PV)