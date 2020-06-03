export
ruleSPSampleListOutNPP

function ruleSPSampleListOutNPP(msg_out::Nothing,msg_s::Message{PointMass,V},msg_w::Message{PointMass, V}) where {V<:VariateType}
    var_type = Univariate
    if length(msg_s.dist.params[:m][1]) > 1
        var_type = Multivariate
    end
    Message(var_type, SampleList, s=deepcopy(msg_s.dist.params[:m]), w=deepcopy(msg_w.dist.params[:m]))
end
