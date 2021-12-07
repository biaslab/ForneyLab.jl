export
ruleSPSampleListOutNPP,
ruleVBSampleListOut

function ruleSPSampleListOutNPP(msg_out::Nothing,
	                            msg_s::Message{PointMass, Multivariate}, # Multivariate, because vectors of samples and weight are passed as parameters
	                            msg_w::Message{PointMass, Multivariate})
    
    s = deepcopy(msg_s.dist.params[:m])
    w = deepcopy(msg_w.dist.params[:m])

    return Message(variateType(s[1]), SampleList, s=s, w=w)
end

function ruleVBSampleListOut(dist_out::Any,
                             dist_s::ProbabilityDistribution{Multivariate, PointMass}, # Multivariate, because vectors of samples and weight are passed as parameters
                             dist_w::ProbabilityDistribution{Multivariate, PointMass})
    
    s = deepcopy(dist_s.params[:m])
    w = deepcopy(dist_w.params[:m])
                         
    return Message(variateType(s[1]), SampleList, s=s, w=w)
end
