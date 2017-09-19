export ruleSPEqualityGaussianMV

ruleSPEqualityGaussianMV(   msg_1::Message{GaussianMeanVariance},
                            msg_2::Message{GaussianMeanVariance},
                            msg_3::Void) =
    Message{GaussianMeanVariance}(msg_1.dist * msg_2.dist)

ruleSPEqualityGaussianMV(   msg_1::Message{GaussianMeanVariance},
                            msg_2::Void,
                            msg_3::Message{GaussianMeanVariance}) =
    Message{GaussianMeanVariance}(msg_1.dist * msg_3.dist)

ruleSPEqualityGaussianMV(   msg_1::Void,
                            msg_2::Message{GaussianMeanVariance},
                            msg_3::Message{GaussianMeanVariance}) =
    Message{GaussianMeanVariance}(msg_2.dist * msg_3.dist)