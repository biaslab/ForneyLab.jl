export
ruleSPNonconjugateInFN,
ruleSPNonconjugateOutNG,
prod!

function ruleSPNonconjugateInFN(msg_out::Message,msg_z::Nothing,g::Function)
    log_grad_g(z) = log(abs(ForwardDiff.derivative(g,z)))
    Message(Univariate, Function, log_pdf=(z)-> log_grad_g(z) + logPdf(msg_out.dist,g(z)))
end

function ruleSPNonconjugateOutNG(msg_out::Nothing,msg_z::Message{F1, V},g::Function) where {F1<:Gaussian, V<:Univariate}
    #z = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_z)
    z_m = msg_z.dist.params[:m]
    z_v = msg_z.dist.params[:v]
    Message(V, Abstract_dist, m=z_m, v=z_v, f=g)
end

const product = Iterators.product
const PIterator = Iterators.ProductIterator

function gaussHermiteQuadrature(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::Array, weights_iter::Array)
   std = sqrt(d.params[:v])
   result = sum((weights_iter./ (2 .^ collect(0:19) .* sqrt(pi))).*g.(d.params[:m] .+ std.*points_iter))
   return result
end

@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, GaussianMeanVariance}=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) where {F2<:Gaussian}

    y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
    y_m = y.params[:m]
    y_v = y.params[:v]
    FB(s) = exp(x.params[:log_pdf](s))
    points_iter, weights_iter = gausshermite(20)
    normalization_constant = gaussHermiteQuadrature(FB,y,points_iter,weights_iter)
    t(s) = s*FB(s)/normalization_constant
    mean = gaussHermiteQuadrature(t,y,points_iter,weights_iter)
    c(s) = FB(s)*(s-mean)^2/normalization_constant
    cov = gaussHermiteQuadrature(c,y,points_iter,weights_iter)
    z.params[:m] = mean
    z.params[:v] = cov

    return z
end

function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, Function},
    z::ProbabilityDistribution{Univariate, Function}=ProbabilityDistribution(Univariate, Function, log_pdf=(z)-> z))

    log_FA = x.params[:log_pdf]
    log_FB = y.params[:log_pdf]

    z.params[:log_pdf] = (s -> log_FA(s) + log_FB(s))

    return z
end

#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::Nonconjugate, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    inbound_messages = String[]
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        if node_interface == entry.interface
            # Ignore inbound message on outbound interface
            push!(inbound_messages, "nothing")
        elseif isa(inbound_interface.node, Clamp)
            # Hard-code outbound message of constant node in schedule
            push!(inbound_messages, messageString(inbound_interface.node))
        else
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbound_messages, "messages[$inbound_idx]")
        end
    end

    push!(inbound_messages, "$(node.g)")

    return inbound_messages
end
