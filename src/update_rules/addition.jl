@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPAdditionOutVGG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing),
                :name          => SPAdditionIn2GGV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{Gaussian}),
                :name          => SPAdditionIn1GVG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPAdditionOutVGP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{Gaussian}),
                :name          => SPAdditionOutVPG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Nothing, Message{Gaussian}),
                :name          => SPAdditionIn1PVG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPAdditionIn1GVP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Message{Gaussian}, Nothing),
                :name          => SPAdditionIn2PGV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Nothing),
                :name          => SPAdditionIn2GPV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPAdditionOutVPP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Nothing),
                :name          => SPAdditionIn2PPV)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPAdditionIn1PVP)
