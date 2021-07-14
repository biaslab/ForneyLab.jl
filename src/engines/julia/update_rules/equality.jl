export
ruleSPEqualityGaussian,
ruleSPEqualityGammaWishart,
ruleSPEqualityBernoulli,
ruleSPEqualityBeta,
ruleSPEqualityCategorical,
ruleSPEqualityDirichlet,
ruleSPEqualityPointMass,
ruleSPEqualityFn,
ruleSPEqualityFnG,
ruleSPEqualityFnFactor,
ruleSPEqualityGFactor,
ruleSPEqualityFactor
# ruleSPEqualityFNFN

ruleSPEqualityGaussian(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Gaussian, F2<:Gaussian} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGaussian(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Gaussian}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGaussian(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Gaussian} = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityGammaWishart(msg_1::Message{F}, msg_2::Message{F}, msg_3::Nothing) where F<:Union{Gamma, Wishart} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGammaWishart(msg_1::Message{F}, msg_2::Nothing, msg_3::Message{F}) where F<:Union{Gamma, Wishart}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGammaWishart(msg_1::Nothing, msg_2::Message{F}, msg_3::Message{F}) where F<:Union{Gamma, Wishart} = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityBernoulli(msg_1::Message{Bernoulli}, msg_2::Message{Bernoulli}, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityBernoulli(msg_1::Message{Bernoulli}, msg_2::Nothing, msg_3::Message{Bernoulli}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityBernoulli(msg_1::Nothing, msg_2::Message{Bernoulli}, msg_3::Message{Bernoulli}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityBeta(msg_1::Message{Beta}, msg_2::Message{Beta}, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityBeta(msg_1::Message{Beta}, msg_2::Nothing, msg_3::Message{Beta}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityBeta(msg_1::Nothing, msg_2::Message{Beta}, msg_3::Message{Beta}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityCategorical(msg_1::Message{Categorical}, msg_2::Message{Categorical}, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityCategorical(msg_1::Message{Categorical}, msg_2::Nothing, msg_3::Message{Categorical}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityCategorical(msg_1::Nothing, msg_2::Message{Categorical}, msg_3::Message{Categorical}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityDirichlet(msg_1::Message{Dirichlet}, msg_2::Message{Dirichlet}, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityDirichlet(msg_1::Message{Dirichlet}, msg_2::Nothing, msg_3::Message{Dirichlet}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityDirichlet(msg_1::Nothing, msg_2::Message{Dirichlet}, msg_3::Message{Dirichlet}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityPointMass(msg_1::Message, msg_2::Message, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityPointMass(msg_1::Message, msg_2::Nothing, msg_3::Message) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityPointMass(msg_1::Nothing, msg_2::Message, msg_3::Message) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityFn(msg_1::Message{Function}, msg_2::Message{Function}, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityFn(msg_1::Message{Function}, msg_2::Nothing, msg_3::Message{Function}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityFn(msg_1::Nothing, msg_2::Message{Function}, msg_3::Message{Function}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityFnG(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Function, F2<:Gaussian} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityFnG(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Function, F2<:Gaussian}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityFnG(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Function, F2<:Gaussian} = Message(prod!(msg_2.dist, msg_3.dist))
ruleSPEqualityFnG(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Gaussian, F2<:Function} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityFnG(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Function}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityFnG(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Function} = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityFnFactor(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Function, F2<:FactorNode} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityFnFactor(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Function, F2<:FactorNode}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityFnFactor(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Function, F2<:FactorNode} = Message(prod!(msg_2.dist, msg_3.dist))
ruleSPEqualityFnFactor(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:FactorNode, F2<:Function} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityFnFactor(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:FactorNode, F2<:Function}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityFnFactor(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:FactorNode, F2<:Function} = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityGFactor(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Gaussian, F2<:FactorNode} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGFactor(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Gaussian, F2<:FactorNode}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGFactor(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Gaussian, F2<:FactorNode} = Message(prod!(msg_2.dist, msg_3.dist))
ruleSPEqualityGFactor(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:FactorNode, F2<:Gaussian} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGFactor(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:FactorNode, F2<:Gaussian}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGFactor(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:FactorNode, F2<:Gaussian} = Message(prod!(msg_2.dist, msg_3.dist))

# # Since CVI sends a FactorNode message, the equality calls ruleSPEqualityGFactor
# # when the other message is Gaussian.
# # Therefore we implement the below rules.
#
# ruleSPEqualityGFactor(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Gaussian, F2<:Gaussian} = Message(prod!(msg_1.dist, msg_2.dist))
# ruleSPEqualityGFactor(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Gaussian}= Message(prod!(msg_1.dist, msg_3.dist))
# ruleSPEqualityGFactor(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Gaussian} = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityFactor(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:FactorNode, F2<:FactorNode} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityFactor(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:FactorNode, F2<:FactorNode}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityFactor(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:FactorNode, F2<:FactorNode} = Message(prod!(msg_2.dist, msg_3.dist))

# ruleSPEqualityFNFN(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:FactorNode, F2<:FactorNode} = Message(prod!(msg_1.dist, msg_2.dist))
# ruleSPEqualityFNFN(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:FactorNode, F2<:FactorNode}= Message(prod!(msg_1.dist, msg_3.dist))
# ruleSPEqualityFNFN(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:FactorNode, F2<:FactorNode} = Message(prod!(msg_2.dist, msg_3.dist))
