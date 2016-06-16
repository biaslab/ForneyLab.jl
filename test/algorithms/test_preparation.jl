#####################
# Unit tests
#####################

facts("Shared preparation methods for inference algorithms") do
    context("injectParameters!() should fill the parameters of a destination distribution with the parameters of a source distribution") do
        source = Gaussian(m=4.0, V=5.0)
        destination = Gaussian()
        ForneyLab.injectParameters!(destination, source)
        @fact destination.m --> 4.0
        @fact destination.V --> 5.0
        @fact is(destination, source) --> false
        @fact_throws ForneyLab.injectParameters!(destination, Delta(4.0))
    end

    context("collectTypeVarNames and resolveTypeVars") do
        # Specify dummy function for tests
        function dummyRule{dims,n_factors}( node::MockNode,
                                    outbound_iface_id,
                                    outbound::Partitioned{MvGaussian{dims},n_factors},
                                    inbound1::Message{MvGaussian{dims}},
                                    inbound2::Message{Partitioned{Gaussian,n_factors}})
            return true
        end

        method = start(methods(dummyRule))

        # collectTypeVarNames should return the names of TypeVar parameters as a Set
        @fact ForneyLab.collectTypeVarNames(method.sig.types[1]) --> Set()
        @fact ForneyLab.collectTypeVarNames(method.sig.types[3]) --> Set([:dims; :n_factors])
        @fact ForneyLab.collectTypeVarNames(method.sig.types[5]) --> Set([:n_factors])

        # resolveTypeVars should resolve the values of TypeVars
        # Build calling signature (vector of argument types) for dummyRule
        # The combination of method and argtypes fully specifies the values of the dims and n_factors parameters of the outbound argument
        argtypes = [Node; Int64; Any; Message{MvGaussian{5}}; Message{Partitioned{Gaussian,3}}]
        @fact ForneyLab.resolveTypeVars(method.sig.types[3], method, argtypes, MockNode()) --> Partitioned{MvGaussian{5},3}
    end

    context("collectAllOutboundTypes() should return outbound types of applicable update rules") do
        FactorGraph()

        call_signature = [TerminalNode, Type{Val{1}}, Any, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, TerminalNode()) --> [Gaussian]

        call_signature = [TerminalNode, Type{Val{1}}, Any, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, TerminalNode([1.0, 1.0])) --> [MvDelta{Float64, 2}]

        call_signature = [EqualityNode, Type{Val{2}}, Any, Message{Gaussian}, Void, Message{Gaussian}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> [Gaussian]

        call_signature = [EqualityNode, Type{Val{3}}, Any, Message{MvGaussian{2}}, Message{MvDelta{Float64, 2}}, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> [MvDelta{Float64, 2}]

        call_signature = [AdditionNode, Type{Val{2}}, Any, Message{MvGaussian{2}}, Void, Message{MvGaussian{2}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, AdditionNode()) --> [MvGaussian{2}]

        call_signature = [AdditionNode, Type{Val{2}}, Any, Message{MvDelta{Float64, 2}}, Void, Message{MvDelta{Float64, 2}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, AdditionNode()) --> [MvDelta{Float64, 2}]

        call_signature = [GaussianNode{Val{:mean}, Val{:precision}}, Type{Val{2}}, Any, MvGaussian{2}, Void, MvGaussian{2}]
        @fact ForneyLab.collectAllOutboundTypes(variationalRule!, call_signature, GaussianNode(form=:precision)) --> [Wishart{2}]

        call_signature = [GaussianNode{Val{:mean}, Val{:precision}}, Type{Val{3}}, Any, NormalGamma, NormalGamma, Void]
        @fact ForneyLab.collectAllOutboundTypes(variationalRule!, call_signature, GaussianNode(form=:precision)) --> [Gaussian]

        call_signature = [GaussianNode{Val{:mean}, Val{:precision}}, Type{Val{1}}, Any, Void, Message{Gamma}, Gaussian]
        @fact ForneyLab.collectAllOutboundTypes(variationalRule!, call_signature, GaussianNode(form=:precision)) --> [StudentsT]

        call_signature = [SigmoidNode, Type{Val{1}}, Any, Message{Gaussian}, Message{Delta{Bool}}]
        @fact ForneyLab.collectAllOutboundTypes(expectationRule!, call_signature, SigmoidNode()) --> [Gaussian]

        # Tests for approximate msg computation rules
        # Fixed outbound distribution type
        call_signature = [EqualityNode, Type{Val{1}}, Any, Any, Message{Gaussian}, Message{StudentsT}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> []
        call_signature = [EqualityNode, Type{Val{1}}, Any, Any, Message{Gaussian}, Message{StudentsT}, Any]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> [Approximation{Gaussian,MomentMatching}]
        # Fixed outbound distribution type and approximation type
        call_signature = [EqualityNode, Type{Val{1}}, Any, Any, Message{Gaussian}, Message{StudentsT}, Type{MomentMatching}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, EqualityNode()) --> [Approximation{Gaussian,MomentMatching}]
    end

    context("collectAllOutboundTypes() should return outbound types for nodes with a fixed gain under inputs of different dimensions") do
        FactorGraph()

        call_signature = [GainNode, Type{Val{1}}, Any, Void, Message{MvGaussian{3}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainNode(gain=rand(3,2))) --> [MvGaussian{2}]

        call_signature = [GainNode, Type{Val{2}}, Any, Message{MvGaussian{2}}, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainNode(gain=rand(3,2))) --> [MvGaussian{3}]

        call_signature = [GainAdditionNode, Type{Val{1}}, Any, Void, Message{MvGaussian{3}}, Message{MvGaussian{3}}]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainAdditionNode(rand(3,2))) --> [MvGaussian{2}]

        call_signature = [GainEqualityNode, Type{Val{3}}, Any, Message{MvGaussian{2}}, Message{MvGaussian{2}}, Void]
        @fact ForneyLab.collectAllOutboundTypes(sumProductRule!, call_signature, GainEqualityNode(rand(3,2))) --> [MvGaussian{3}]
    end

    FactorGraph()

    context("inferOutboundType!() should infer the correct outbound type for a terminal node") do
        entry = ScheduleEntry(TerminalNode(Gaussian()), 1, sumProductRule!)
        entry.inbound_types = [Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gaussian
    end

    context("inferOutboundType!() should infer the correct outbound type for a node") do
        entry = ScheduleEntry(AdditionNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{Gaussian}, Message{Gaussian}, Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gaussian
    end

    context("inferOutboundType!() should infer the correct outbound type for a Gaussian node with vmp") do
        entry = ScheduleEntry(GaussianNode(form=:precision), 2, variationalRule!)
        entry.inbound_types = [Gaussian, Void, Gaussian]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gamma
    end

    context("inferOutboundType!() should consider approximate rules") do
        # Just one approximate rule available
        entry = ScheduleEntry(EqualityNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{Gaussian}, Message{StudentsT}, Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gaussian
        @fact entry.approximation --> MomentMatching

        # Just one approximate rule available, fixed outbound type
        entry = ScheduleEntry(EqualityNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{Gaussian}, Message{StudentsT}, Void]
        entry.outbound_type = Gaussian
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gaussian
        @fact entry.approximation --> MomentMatching

        # Multiple approximate rules available, fixed outbound type
        entry = ScheduleEntry(ExponentialNode(), 2, sumProductRule!)
        entry.inbound_types = [Message{Gaussian}, Void]
        entry.outbound_type = Gamma
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gamma
        @fact isdefined(entry, :approximation) --> true # we don't care which approximation has been chosen


        # Multiple approximate rules available, fixed outbound type and approximation type
        entry = ScheduleEntry(ExponentialNode(), 2, sumProductRule!)
        entry.inbound_types = [Message{Gaussian}, Void]
        entry.outbound_type = Gamma
        entry.approximation = MomentMatching
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Gamma
        @fact entry.approximation --> MomentMatching
    end

    context("inferOutboundType!() should consider one-by-one processing of Partitioned factors") do
        # Only Partitioned inbounds
        entry = ScheduleEntry(EqualityNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{Partitioned{Gaussian,3}}, Message{Partitioned{Gaussian,3}}, Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Partitioned{Gaussian,3}

        # Mismatch in number of factors
        entry = ScheduleEntry(EqualityNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{Partitioned{Gaussian,2}}, Message{Partitioned{Gaussian,3}}, Void]
        @fact_throws ForneyLab.inferOutboundType!(entry)

        # Partitioned inbound mixed with regular inbound
        entry = ScheduleEntry(EqualityNode(), 3, sumProductRule!)
        entry.inbound_types = [Message{Gaussian}, Message{Partitioned{Gaussian,3}}, Void]
        ForneyLab.inferOutboundType!(entry)
        @fact entry.outbound_type --> Partitioned{Gaussian,3}
    end

    context("buildExecute!() should pre-compile the execute field of the schedule entry") do
        node = AdditionNode()
        node.i[:out].message = Message(Gaussian())
        entry = ScheduleEntry(node, 3, sumProductRule!)
        @fact isdefined(entry, :execute) --> false
        ForneyLab.buildExecute!(entry, [Message(Gaussian()), Message(Gaussian()), nothing])
        @fact typeof(entry.execute) --> Function
    end
end
