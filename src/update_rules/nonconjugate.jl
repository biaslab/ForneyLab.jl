@sumProductRule(:node_type     => Nonconjugate,
                :outbound_type => Message{Function},
                :inbound_types => (Message{Union{Bernoulli, Beta, Categorical, Dirichlet, Gaussian, Gamma, LogNormal, Poisson, Wishart}}, Nothing),
                :name          => SPNonconjugateInFN)

@sumProductRule(:node_type     => Nonconjugate,
                :outbound_type => Message{Abstract_dist},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonconjugateOutNG)
