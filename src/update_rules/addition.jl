@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{Gaussian}, Message{Gaussian}),
                :name          => SPAdditionOutVGG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Void),
                :name          => SPAdditionIn2GGV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Void, Message{Gaussian}),
                :name          => SPAdditionIn1GVG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPAdditionOutVGP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPAdditionOutVPG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Void, Message{Gaussian}),
                :name          => SPAdditionIn1PVG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPAdditionIn1GVP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Message{Gaussian}, Void),
                :name          => SPAdditionIn2PGV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Void),
                :name          => SPAdditionIn2GPV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPAdditionOutVPP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Void),
                :name          => SPAdditionIn2PPV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Void, Message{PointMass}),
                :name          => SPAdditionIn1PVP)
