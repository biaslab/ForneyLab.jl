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

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{SampleList}, Message{PointMass}),
                :name          => SPAdditionOutNSP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{PointMass}, Message{SampleList}),
                :name          => SPAdditionOutNPS)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{PointMass}, Nothing, Message{SampleList}),
                :name          => SPAdditionIn1PNS)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{SampleList}, Nothing, Message{PointMass}),
                :name          => SPAdditionIn1SNP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{PointMass}, Message{SampleList}, Nothing),
                :name          => SPAdditionIn2PSN)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{SampleList}, Message{PointMass}, Nothing),
                :name          => SPAdditionIn2SPN)
                
@marginalRule(:node_type     => Addition,
              :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}), # "Nothing" indicates marginalization over outbound edge
              :name          => MAdditionNGG)

