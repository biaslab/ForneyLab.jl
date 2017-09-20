@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Void),
                :name          => SPAdditionGGV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{Gaussian}),
                :name          => SPAdditionGVG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{Gaussian}),
                :name          => SPAdditionVGG)

# TODO: add other combinations
@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPAdditionVGP)