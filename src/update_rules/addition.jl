@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{Gaussian}),
                :name          => SPAdditionOutGG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Void),
                :name          => SPAdditionIn2GG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{Gaussian}),
                :name          => SPAdditionIn1GG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPAdditionOutGP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPAdditionOutPG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Void, Message{Gaussian}),
                :name          => SPAdditionIn1PG)