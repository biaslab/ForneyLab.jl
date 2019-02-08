export
ruleSPGainEqualityOutVGGP,
ruleSPGainEqualityIn1GVGP

function ruleSPGainEqualityOutVGGP(msg_out::Nothing,
                                   msg_in1::Message{F, Multivariate},
                                   msg_in2::Message{F, V1},
                                   gain::Message{PointMass, V2}) where {F<:Gaussian, V1, V2}

   isa(msg_in2, Message{Gaussian, Univariate}) ? gainfactor = gain.dist.params[:m]' : gainfactor = gain.dist.params[:m]
   in1_m = unsafeMean(msg_in1.dist); in2_m = unsafeMean(msg_in2.dist)
   in1_w = unsafePrecision(msg_in1.dist); in2_w = unsafePrecision(msg_in2.dist)
   W_out = in1_w + gainfactor*in2_w*gainfactor'
   out_m = (W_out)^-1*(in1_w*in1_m + gainfactor*in2_w*in2_m)
   return Message(Multivariate, GaussianMeanPrecision, m=out_m, w=W_out)

end

function ruleSPGainEqualityIn1GVGP(msg_out::Message{F, Multivariate},
                                   msg_in1::Nothing,
                                   msg_in2::Message{F, V1},
                                   gain::Message{PointMass, V2}) where {F<:Gaussian, V1, V2}

   return ruleSPGainEqualityOutVGGP(msg_in1, msg_out, msg_in2, gain)
end
