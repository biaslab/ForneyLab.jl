@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Univariate, Void, Univariate, Univariate, Univariate, Univariate),
                 :name          => VBGaussianMixtureM1)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Univariate{Gamma}},
                 :inbound_types => (Univariate, Univariate, Void, Univariate, Univariate, Univariate),
                 :name          => VBGaussianMixtureW1)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Univariate, Univariate, Univariate, Void, Univariate, Univariate),
                 :name          => VBGaussianMixtureM2)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Univariate{Gamma}},
                 :inbound_types => (Univariate, Univariate, Univariate, Univariate, Void, Univariate),
                 :name          => VBGaussianMixtureW2)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Univariate{Bernoulli}},
                 :inbound_types => (Univariate, Univariate, Univariate, Univariate, Univariate, Void),
                 :name          => VBGaussianMixtureZ)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Void, Univariate, Univariate, Univariate, Univariate, Univariate),
                 :name          => VBGaussianMixtureOut)