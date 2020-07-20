export
ruleSPSampleListOutNPP

function ruleSPSampleListOutNPP(msg_out::Nothing,
	                            msg_s::Message{PointMass, Multivariate}, # Multivariate, because vectors of samples and weight are passed as parameters
	                            msg_w::Message{PointMass, Multivariate})
    
    # Extract the variate type from the first sample
    s = msg_s.dist.params[:m][1]
    if isa(s, Number)
	    V = Univariate
	elseif isa(s, AbstractVector)
        V = Multivariate
    elseif isa(s, AbstractMatrix)
       	V = MatrixVariate
    else
    	error("Unexpected sample type: $(typeof(s))")
    end

    return Message(V, SampleList, s=deepcopy(msg_s.dist.params[:m]), w=deepcopy(msg_w.dist.params[:m]))
end
