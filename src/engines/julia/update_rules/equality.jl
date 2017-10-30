export
ruleSPEqualityGaussian, 
ruleSPEqualityGamma

function ruleSPEqualityGaussian(msg_1::Message{Univariate{Gaussian}},
                                msg_2::Message{Univariate{Gaussian}},
                                msg_3::Void)
    
    outbound_dist = Univariate(Gaussian, xi=0.0, w=1.0)
    prod!(msg_1.dist, msg_2.dist, outbound_dist)

    return Message{Univariate{Gaussian}}(outbound_dist)
end
ruleSPEqualityGaussian(msg_1::Message{Univariate{Gaussian}}, msg_2::Void, msg_3::Message{Univariate{Gaussian}}) = ruleSPEqualityGaussian(msg_1, msg_3, msg_2)
ruleSPEqualityGaussian(msg_1::Void, msg_2::Message{Univariate{Gaussian}}, msg_3::Message{Univariate{Gaussian}}) = ruleSPEqualityGaussian(msg_2, msg_3, msg_1)

function ruleSPEqualityGamma(   msg_1::Message{Univariate{Gamma}},
                                msg_2::Message{Univariate{Gamma}},
                                msg_3::Void)
    
    outbound_dist = Univariate(Gamma, a=0.0, b=0.0)
    prod!(msg_1.dist, msg_2.dist, outbound_dist)

    return Message{Univariate{Gamma}}(outbound_dist)
end
ruleSPEqualityGamma(msg_1::Message{Univariate{Gamma}}, msg_2::Void, msg_3::Message{Univariate{Gamma}}) = ruleSPEqualityGamma(msg_1, msg_3, msg_2)
ruleSPEqualityGamma(msg_1::Void, msg_2::Message{Univariate{Gamma}}, msg_3::Message{Univariate{Gamma}}) = ruleSPEqualityGamma(msg_2, msg_3, msg_1)