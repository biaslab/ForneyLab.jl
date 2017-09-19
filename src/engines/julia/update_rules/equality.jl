export ruleSPEqualityGaussian

function ruleSPEqualityGaussian{T<:Gaussian, U<:Gaussian}(  msg_1::Message{T},
                                                            msg_2::Message{U},
                                                            msg_3::Void)
    
    output_dist = ProbabilityDistribution(GaussianWeightedMeanPrecision)
    prod!(msg_1.dist, msg_2.dist, output_dist)

    return Message{GaussianWeightedMeanPrecision}(output_dist)
end

ruleSPEqualityGaussian{T<:Gaussian, U<:Gaussian}(msg_1::Message{T}, msg_2::Void, msg_3::Message{U}) = ruleSPEqualityGaussian(msg_1, msg_3, msg_2)
ruleSPEqualityGaussian{T<:Gaussian, U<:Gaussian}(msg_1::Void, msg_2::Message{T}, msg_3::Message{U}) = ruleSPEqualityGaussian(msg_2, msg_3, msg_1)