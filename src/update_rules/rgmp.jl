@sumProductRule(:node_type     => RGMP,
                :outbound_type => Message{Function},
                :inbound_types => (Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing),
                :name          => SPRGMPInFN)

@sumProductRule(:node_type     => RGMP,
                :outbound_type => Message{RGMP_dist},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPRGMPOutNG)
