@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Void, Message{Univariate{Gaussian}}, Message{Univariate{PointMass}}),
                :name          => SPMultiplicationOutGP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Message{Univariate{Gaussian}}, Void, Message{Univariate{PointMass}}),
                :name          => SPMultiplicationInGP)