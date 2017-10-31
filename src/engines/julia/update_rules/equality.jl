export
ruleSPEqualityGaussian, 
ruleSPEqualityGamma

function ruleSPEqualityGaussian(msg_1::Message{Gaussian, Univariate},
                                msg_2::Message{Gaussian, Univariate},
                                msg_3::Void)
    
    outbound_dist = Univariate(Gaussian, xi=0.0, w=1.0)
    prod!(msg_1.dist, msg_2.dist, outbound_dist)

    return Message(outbound_dist)
end
ruleSPEqualityGaussian(msg_1::Message{Gaussian, Univariate}, msg_2::Void, msg_3::Message{Gaussian, Univariate}) = ruleSPEqualityGaussian(msg_1, msg_3, msg_2)
ruleSPEqualityGaussian(msg_1::Void, msg_2::Message{Gaussian, Univariate}, msg_3::Message{Gaussian, Univariate}) = ruleSPEqualityGaussian(msg_2, msg_3, msg_1)

function ruleSPEqualityGamma(   msg_1::Message{Gamma, Univariate},
                                msg_2::Message{Gamma, Univariate},
                                msg_3::Void)
    
    outbound_dist = Univariate(Gamma, a=0.0, b=0.0)
    prod!(msg_1.dist, msg_2.dist, outbound_dist)

    return Message(outbound_dist)
end
ruleSPEqualityGamma(msg_1::Message{Gamma, Univariate}, msg_2::Void, msg_3::Message{Gamma, Univariate}) = ruleSPEqualityGamma(msg_1, msg_3, msg_2)
ruleSPEqualityGamma(msg_1::Void, msg_2::Message{Gamma, Univariate}, msg_3::Message{Gamma, Univariate}) = ruleSPEqualityGamma(msg_2, msg_3, msg_1)