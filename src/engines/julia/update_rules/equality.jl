export
ruleSPEqualityGaussian,
ruleSPEqualityGammaWishart,
ruleSPEqualityBernoulli,
ruleSPEqualityBeta,
ruleSPEqualityCategorical,
ruleSPEqualityDirichlet,
ruleSPEqualityPointMass,
ruleSPEqualityRGMP,
ruleSPEqualityGaussianRGMP

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

ruleSPEqualityRGMP(msg_1::Message{Function}, msg_2::Message{Function}, msg_3::Nothing) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityRGMP(msg_1::Message{Function}, msg_2::Nothing, msg_3::Message{Function}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityRGMP(msg_1::Nothing, msg_2::Message{Function}, msg_3::Message{Function}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityGaussianRGMP(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Function, F2<:Gaussian} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGaussianRGMP(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Function, F2<:Gaussian}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGaussianRGMP(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Function, F2<:Gaussian} = Message(prod!(msg_2.dist, msg_3.dist))
ruleSPEqualityGaussianRGMP(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Nothing) where {F1<:Gaussian, F2<:Function} = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGaussianRGMP(msg_1::Message{F1}, msg_2::Nothing, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Function}= Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGaussianRGMP(msg_1::Nothing, msg_2::Message{F1}, msg_3::Message{F2}) where {F1<:Gaussian, F2<:Function} = Message(prod!(msg_2.dist, msg_3.dist))
