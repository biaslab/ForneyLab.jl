export
ruleSPUnitEqualityOutVGG,
ruleSPUnitEqualityIn1GVG

function ruleSPUnitEqualityOutVGG(msg_out::Nothing,
                                  msg_in1::Message{F, Multivariate},
                                  msg_in2::Message{F, Univariate}) where F<:Gaussian

   in1_m = unsafeMean(msg_in1.dist)
   in2_m = unsafeMean(msg_in2.dist)
   dim = length(in1_m)
   uvector = zeros(dim); uvector[1] = 1;
   W_out = msg_in1.dist.params[:w] + uvector*msg_in2.dist.params[:w]*uvector'
   out_m = (W_out)^-1*(msg_in1.dist.params[:w]*msg_in1.dist.params[:m] + uvector*msg_in2.dist.params[:w]*msg_in2.dist.params[:m])
   return Message(Multivariate, GaussianMeanPrecision, m=out_m, w=W_out)

end

function ruleSPUnitEqualityIn1GVG(msg_out::Message{F, Multivariate},
                                  msg_in1::Nothing,
                                  msg_in2::Message{F, Univariate}) where F<:Gaussian

   return ruleSPUnitEqualityOutVGG(msg_in1, msg_out, msg_in2)
end
