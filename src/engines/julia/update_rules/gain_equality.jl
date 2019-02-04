export
ruleVBGammaOut,
ruleSPGammaOutVPP

ruleSPGainEqualityOutVGG(msg_out::Nothing,
                         msg_in1::Message{PointMass, Univariate},
                         msg_in2::Message{PointMass, Univariate}) =
                         Message(Multivariate, Gaussian, a=deepcopy(msg_in1.dist.params[:m]), b=deepcopy(msg_in2.dist.params[:m]))



function ruleSPGainEqualityOutVGG(msg_out::Nothing,
                                  msg_in1::Message{F, Multivariate},
                                  msg_in2::Message{F, Univariate}) where F<:Gaussian

   in1_m = unsafeMean(msg_in1)
   in2_m = unsafeMean(msg_in2)
   dim = length(in1_m)
   W_out = msg_in1.dist.params[:w] + uvector(dim)*msg_in1.dist.params[:w]*uvector(dim)'
   out_m = (W_out)^-1*(msg_in1.dist.params[:w]*msg_in1.dist.params[:xi] + uvector(dim)*msg_in2.dist.params[:w]*msg_in2.dist.params[:w])

   return Message(Univariate, GaussianMeanPrecision, m=out_m, w=W_out)

end
