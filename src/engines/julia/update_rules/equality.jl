export
ruleSPEqualityGaussian, 
ruleSPEqualityGammaWishart,
ruleSPEqualityBernoulli,
ruleSPEqualityBeta,
ruleSPEqualityCategorical,
ruleSPEqualityDirichlet,
ruleSPEqualityPointMass

ruleSPEqualityGaussian{F1<:Gaussian, F2<:Gaussian}(msg_1::Message{F1}, msg_2::Message{F2}, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGaussian{F1<:Gaussian, F2<:Gaussian}(msg_1::Message{F1}, msg_2::Void, msg_3::Message{F2}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGaussian{F1<:Gaussian, F2<:Gaussian}(msg_1::Void, msg_2::Message{F1}, msg_3::Message{F2}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityGammaWishart{F<:Union{Gamma, Wishart}}(msg_1::Message{F}, msg_2::Message{F}, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityGammaWishart{F<:Union{Gamma, Wishart}}(msg_1::Message{F}, msg_2::Void, msg_3::Message{F}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityGammaWishart{F<:Union{Gamma, Wishart}}(msg_1::Void, msg_2::Message{F}, msg_3::Message{F}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityBernoulli(msg_1::Message{Bernoulli}, msg_2::Message{Bernoulli}, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityBernoulli(msg_1::Message{Bernoulli}, msg_2::Void, msg_3::Message{Bernoulli}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityBernoulli(msg_1::Void, msg_2::Message{Bernoulli}, msg_3::Message{Bernoulli}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityBeta(msg_1::Message{Beta}, msg_2::Message{Beta}, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityBeta(msg_1::Message{Beta}, msg_2::Void, msg_3::Message{Beta}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityBeta(msg_1::Void, msg_2::Message{Beta}, msg_3::Message{Beta}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityCategorical(msg_1::Message{Categorical}, msg_2::Message{Categorical}, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityCategorical(msg_1::Message{Categorical}, msg_2::Void, msg_3::Message{Categorical}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityCategorical(msg_1::Void, msg_2::Message{Categorical}, msg_3::Message{Categorical}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityDirichlet(msg_1::Message{Dirichlet}, msg_2::Message{Dirichlet}, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityDirichlet(msg_1::Message{Dirichlet}, msg_2::Void, msg_3::Message{Dirichlet}) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityDirichlet(msg_1::Void, msg_2::Message{Dirichlet}, msg_3::Message{Dirichlet}) = Message(prod!(msg_2.dist, msg_3.dist))

ruleSPEqualityPointMass(msg_1::Message, msg_2::Message, msg_3::Void) = Message(prod!(msg_1.dist, msg_2.dist))
ruleSPEqualityPointMass(msg_1::Message, msg_2::Void, msg_3::Message) = Message(prod!(msg_1.dist, msg_3.dist))
ruleSPEqualityPointMass(msg_1::Void, msg_2::Message, msg_3::Message) = Message(prod!(msg_2.dist, msg_3.dist))