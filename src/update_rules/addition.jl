@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Void),
                :name          => SPAdditionGGV)