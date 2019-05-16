@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPAdditionOutNGG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing),
                :name          => SPAdditionIn2GGN)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{Gaussian}),
                :name          => SPAdditionIn1GNG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPAdditionOutNGP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{Gaussian}),
                :name          => SPAdditionOutNPG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Nothing, Message{Gaussian}),
                :name          => SPAdditionIn1PNG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPAdditionIn1GNP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Message{Gaussian}, Nothing),
                :name          => SPAdditionIn2PGN)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Nothing),
                :name          => SPAdditionIn2GPN)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPAdditionOutNPP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Nothing),
                :name          => SPAdditionIn2PPN)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPAdditionIn1PNP)
