@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{LogNormal},
                :inbound_types => (Void, Message{Gaussian}),
                :name          => SPExponentialOutVG)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{LogNormal}, Void),
                :name          => SPExponentialIn1LV)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{PointMass},
                :inbound_types => (Void, Message{PointMass}),
                :name          => SPExponentialOutVP)

@sumProductRule(:node_type     => Exponential,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Void),
                :name          => SPExponentialIn1PV)