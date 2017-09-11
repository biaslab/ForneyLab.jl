@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Void),
                :name          => SPGaussianMeanVariancePPV)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceVPP)

# TODO: this function will not work for a GaussianMeanPrecision input on interface 1
@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{GaussianMeanVariance}, Message{PointMass}, Void),
                :name          => SPGaussianMeanVarianceGPV)

