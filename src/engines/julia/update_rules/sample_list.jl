export
ruleSPSampleListOutNPP

function ruleSPSampleListOutNPP(msg_out::Nothing,
	                            msg_s::Message{PointMass, V},
	                            msg_w::Message{PointMass, V}) where V<:VariateType
    
    d = length(msg_s.dist.params[:m][1])
    if d == 1
	    var_type = Univariate
	else
        var_type = Multivariate
    end

    return Message(var_type, SampleList, s=deepcopy(msg_s.dist.params[:m]), w=deepcopy(msg_w.dist.params[:m]))
end
