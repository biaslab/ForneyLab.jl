export
ruleSPRGMPInFN,
ruleSPRGMPOutNG,
prod!

function ruleSPRGMPInFN(msg_out::Message,msg_z::Nothing,g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64)
    log_grad_g(z) = log(abs(ForwardDiff.derivative(g,z)))
    Message(Univariate, Function, log_pdf=(z)-> log_grad_g(z) + logPdf(msg_out.dist,g(z)), eta_a=lr, eta_b=lr, num_epoch=num_epochs1, num_epoch2=num_epochs2)
end

function ruleSPRGMPOutNG(msg_out::Nothing,msg_z::Message{F1, V},g::Function, num_epochs1::Int64, num_epochs2::Int64, lr::Float64) where {F1<:Gaussian, V<:Univariate}
    #z = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_z)
    z_m = msg_z.dist.params[:m]
    z_v = msg_z.dist.params[:v]
    Message(V, RGMP_dist, m=z_m, v=z_v, f=g)
end

const product = Iterators.product
const PIterator = Iterators.ProductIterator

function generate_multidim_points(n::Int, p::Int)
   sigma_points, sigma_weights = gausshermite(p)
   points_iter = product(repeat([sigma_points],n)...)
   weights_iter = product(repeat([sigma_weights],n)...)
   return points_iter, weights_iter
end
function gaussHermiteCubature(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::PIterator, weights_iter::PIterator)
   result = 0.0
   # println(d.params[:v])
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result += weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
function gaussHermiteCubature1D(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::Array, weights_iter::Array)
   result = 0.0
   # println(d.params[:v])
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result += weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
function gaussHermiteCubatureMean(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::PIterator, weights_iter::PIterator)
   result = zeros(dims(d))
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result = result + weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
function gaussHermiteCubatureCov(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::PIterator, weights_iter::PIterator)
   result = zeros(dims(d),dims(d))
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result = result + weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end


@symmetrical function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, F2},
    z::ProbabilityDistribution{Univariate, GaussianMeanVariance}=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) where {F2<:Gaussian}

    y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
    y_m = y.params[:m]
    y_v = y.params[:v]
    log_FB(s) = exp(x.params[:log_pdf](s))
    dim = dims(y)
    #points_iter, weights_iter = generate_multidim_points(dim, 10)
    points_iter, weights_iter = gausshermite(100)
    normalization_constant = gaussHermiteCubature1D(log_FB,y,points_iter,weights_iter)
    t(s) = s.*log_FB(s)./normalization_constant
    #mean = gaussHermiteCubatureMean(t,y,points_iter,weights_iter)
    mean = gaussHermiteCubature1D(t,y,points_iter,weights_iter)
    c(s) = log_FB(s) .* (s-mean)*transpose(s-mean)./normalization_constant
    #cov = gaussHermiteCubatureCov(d,y,points_iter,weights_iter)
    cov = gaussHermiteCubature1D(c,y,points_iter,weights_iter)
    z.params[:m] = mean
    z.params[:v] = cov

    return z
end

function prod!(
    x::ProbabilityDistribution{Univariate, Function},
    y::ProbabilityDistribution{Univariate, Function},
    z::ProbabilityDistribution{Univariate, Function}=ProbabilityDistribution(Univariate, Function, log_pdf=(z)-> z, eta_a=0.01, eta_b=0.01, num_epoch=1, num_epoch2=1000))

    log_FA = x.params[:log_pdf]
    log_FB = y.params[:log_pdf]
    #approximate the posterior with Gaussian distribution family
    eta_a = x.params[:eta_a]
    eta_b = x.params[:eta_b]
    num_epoch = x.params[:num_epoch]
    num_epoch2 = x.params[:num_epoch2]

    z.params[:log_pdf] = (s -> log_FA(s) + log_FB(s))
    z.params[:eta_a] = eta_a
    z.params[:eta_b] = eta_b
    z.params[:num_epoch] = num_epoch
    z.params[:num_epoch2] = num_epoch2

    return z
end

#--------------------------
# Custom inbounds collector
#--------------------------

function collectSumProductNodeInbounds(node::RGMP, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
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
    push!(inbound_messages, "$(node.num_epochs1)")
    push!(inbound_messages, "$(node.num_epochs2)")
    push!(inbound_messages, "$(node.lr)")

    return inbound_messages
end
