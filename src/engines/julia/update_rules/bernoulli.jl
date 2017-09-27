export ruleSPBernoulliPV

ruleSPBernoulliPV(msg_in::Message{PointMass}, msg_out::Void) = Message(Bernoulli, p=msg_in.dist.params[:m])
