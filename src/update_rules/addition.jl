@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{Gaussian}),
                :name          => SPAdditionOutVGG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Void),
                :name          => SPAdditionIn2GGV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{Gaussian}),
                :name          => SPAdditionIn1GVG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPAdditionOutVGP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPAdditionOutVPG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Void, Message{Gaussian}),
                :name          => SPAdditionIn1PVG)