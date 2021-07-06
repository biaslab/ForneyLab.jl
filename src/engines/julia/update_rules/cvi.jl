export ruleSPCVIOutVD, ruleSPCVIIn1MV, ruleSPCVIOutVDX, ruleSPCVIInX

function ruleSPCVIIn1MV(node_id::Symbol,
                        msg_out::Message{<:FactorFunction, <:VariateType},
                        msg_in::Message{<:FactorNode, <:VariateType})

    @show msg_in
    @show typeof(msg_in)
    @show msg_out
    @show node_id
    msg_in
end

function ruleSPCVIIn1MV(node_id::Symbol,
                        msg_out::Message{<:FactorFunction, <:VariateType},
                        msg_in::Message{<:Gaussian, Univariate})

    thenode = currentGraph().nodes[node_id]

    η = deepcopy(naturalParams(msg_in.dist))
    λ = deepcopy(η)

    logp_nc(z) = logPdf(msg_out.dist, thenode.g(z))
    df_m(z) = ForwardDiff.derivative(logp_nc,z)
    df_v(z) = 0.5*ForwardDiff.derivative(df_m,z)
    for i=1:thenode.num_iterations
        q = standardDist(msg_in.dist,λ)
        z_s = sample(q)
        df_μ1 = df_m(z_s) - 2*df_v(z_s)*mean(q)
        df_μ2 = df_v(z_s)
        ∇f = [df_μ1, df_μ2]
        ∇ = λ - η - ∇f
        update!(thenode.opt,λ,∇)
    end

    return standardMessage(msg_in.dist,λ-η)

end

function ruleSPCVIIn1MV(node_id::Symbol,
                        #msg_out::Message{GaussianWeightedMeanPrecision, Multivariate},
                        msg_out::Message{<:FactorFunction, <:VariateType},
                        msg_in::Message{<:Gaussian, Multivariate})

    thenode = currentGraph().nodes[node_id]

    η = deepcopy(naturalParams(msg_in.dist))
    λ = deepcopy(η)

    logp_nc(z) = logPdf(msg_out.dist, thenode.g(z))
    df_m(z) = ForwardDiff.gradient(logp_nc,z)
    df_v(z) = 0.5*ForwardDiff.jacobian(df_m,z)
    for i=1:thenode.num_iterations
        q = standardDist(msg_in.dist,λ)
        z_s = sample(q)
        df_μ1 = df_m(z_s) - 2*df_v(z_s)*mean(q)
        df_μ2 = df_v(z_s)
        ∇f = [df_μ1; vec(df_μ2)]
        λ_old = deepcopy(λ)
        ∇ = λ - η - ∇f
        update!(thenode.opt,λ,∇)
        if isProper(standardDist(msg_in.dist,λ)) == false
            @show node_id
            @show i
            λ = λ_old
        end
    end

    return standardMessage(msg_in.dist,λ-η)
end

function ruleSPCVIOutVD(node_id::Symbol,
                        msg_out::Any,
                        msg_in::ProbabilityDistribution)

    thenode = currentGraph().nodes[node_id]

    samples = thenode.g.(sample(msg_in, thenode.num_samples))
    weights = ones(thenode.num_samples)/thenode.num_samples

    if length(samples[1]) == 1
        variate = Univariate
    else
        variate = Multivariate
    end

    q=ProbabilityDistribution(variate, SampleList, s=samples, w=weights)
    q.params[:entropy] = 0

    return Message(variate,SetSampleList,q=q,node_id=node_id)

end

function ruleSPCVIOutVD(node_id::Symbol,
                        msg_out::Any,
                        msg_in::Message)

    thenode = currentGraph().nodes[node_id]

    samples = thenode.g.(sample(msg_in.dist, thenode.num_samples))
    weights = ones(thenode.num_samples)/thenode.num_samples

    if length(samples[1]) == 1
        variate = Univariate
    else
        variate = Multivariate
    end

    q=ProbabilityDistribution(variate, SampleList, s=samples, w=weights)
    q.params[:entropy] = 0

    return Message(variate,SetSampleList,q=q,node_id=node_id)

end

function ruleSPCVIInX(node_id::Symbol,
                      inx::Int64,
                      msg_out::Message{<:FactorFunction, <:VariateType},
                      msgs_in::Vararg{Union{Message,ProbabilityDistribution}})

    @show msgs_in
    msgs_in[inx]
end

function ruleSPCVIOutVDX(node_id::Symbol,
                         msg_out::Any,
                         msgs_in::Vararg{ProbabilityDistribution})

    thenode = currentGraph().nodes[node_id]

    samples_in = [sample(msg_in, thenode.num_samples) for msg_in in msgs_in]

    samples = thenode.g.(samples_in...)
    weights = ones(thenode.num_samples)/thenode.num_samples

    if length(samples[1]) == 1
        variate = Univariate
    else
        variate = Multivariate
    end

    q=ProbabilityDistribution(variate, SampleList, s=samples, w=weights)
    q.params[:entropy] = 0

    return Message(variate,SetSampleList,q=q,node_id=node_id)

end


#---------------------------
# Custom inbounds collectors
#---------------------------

function collectSumProductNodeInbounds(node::CVI, entry::ScheduleEntry)
    interface_to_schedule_entry = current_inference_algorithm.interface_to_schedule_entry
    target_to_marginal_entry = current_inference_algorithm.target_to_marginal_entry

    inbounds = Any[]

    push!(inbounds, node.id)

    multi_in = (length(node.interfaces) > 2) # Boolean to indicate a multi-inbound nonlinear node
    inx = findfirst(isequal(entry.interface), node.interfaces) - 1 # Find number of inbound interface; 0 for outbound

    if (inx > 0) && multi_in # Multi-inbound backward rule
        push!(inbounds, Dict{Symbol, Any}(:inx => inx, # Push inbound identifier
                                          :keyword => false))
    end

    for node_interface in node.interfaces
        inbound_interface = ultimatePartner(node_interface)

        if (node_interface == entry.interface != node.interfaces[1])
            # Collect the incoming message

            if typeof(inbound_interface.node) == Equality
                # If CVI is connected to Equality node, incoming message is often not available
                # in time series models. The rules here allow us to use the message coming to
                # interface 1 of equality node.
                if inbound_interface == inbound_interface.node.interfaces[2]
                    push!(inbounds, interface_to_schedule_entry[ultimatePartner(inbound_interface.node.interfaces[1])])
                else
                    push!(inbounds, interface_to_schedule_entry[inbound_interface])
                end
            else
                push!(inbounds, interface_to_schedule_entry[inbound_interface])
            end
            #push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif (node_interface == node.interfaces[1] != entry.interface)
            # Collect the BP message from out interface
            push!(inbounds, interface_to_schedule_entry[inbound_interface])
        elseif node_interface === entry.interface
            # Ignore marginal of outbound edge
            push!(inbounds, nothing)
        elseif isClamped(inbound_interface)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, assembleClamp!(inbound_interface.node, ProbabilityDistribution))
        else
            # Collect entry from marginal schedule
            try
                push!(inbounds, target_to_marginal_entry[node_interface.edge.variable])
            catch
                # This rule is useful for the last time step in a time series model with Structured VMP
                push!(inbounds, interface_to_schedule_entry[inbound_interface])
            end
        end
    end

    return inbounds
end
