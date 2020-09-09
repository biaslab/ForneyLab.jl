export
ruleSPNonlinearEOutNG,
ruleSPNonlinearEOutNGX,
ruleSPNonlinearEIn1GG,
ruleSPNonlinearEInGX,
ruleMNonlinearEInGX

"""
Concatenate a vector of vectors and return with original dimensions (for splitting)
"""
function concatenate(xs::Vector{Vector{Float64}})
    ds = [length(x_k) for x_k in xs] # Extract dimensions
    x = vcat(xs...)

    return (x, ds)
end

"""
Return local linearization of g around expansion point x_hat
"""
function localLinearization(V::Type{Univariate}, g::Function, x_hat::Float64)
    a = ForwardDiff.derivative(g, x_hat)
    b = g(x_hat) - a*x_hat

    return (a, b)
end

function localLinearization(V::Type{Multivariate}, g::Function, x_hat::Vector{Float64})
    A = ForwardDiff.jacobian(g, x_hat)
    b = g(x_hat) - A*x_hat

    return (A, b)
end

function localLinearization(V::Type{Univariate}, g::Function, x_hat::Vector{Float64})
    g_unpacked(x::Vector) = g(x...)
    A = ForwardDiff.gradient(g_unpacked, x_hat)'
    b = g(x_hat...) - A*x_hat

    return (A, b)
end

function localLinearization(V::Type{Multivariate}, g::Function, x_hat::Vector{Vector{Float64}})
    (x_cat, ds) = concatenate(x_hat)
    g_unpacked(x::Vector) = g(split(x, ds)...)
    A = ForwardDiff.jacobian(g_unpacked, x_cat)
    b = g(x_hat...) - A*x_cat

    return (A, b)
end


#-----------------------
# Extended Update Rules
#-----------------------

# Forward rule
function ruleSPNonlinearEOutNG(g::Function,
                               msg_out::Nothing,
                               msg_in1::Message{F, V}) where {F<:Gaussian, V<:VariateType}
    
    (m_in1, V_in1) = unsafeMeanCov(msg_in1.dist)
    (A, b) = localLinearization(V, g, m_in1)

    return Message(V, GaussianMeanVariance, m=A*m_in1 + b, v=A*V_in1*A')
end

# Multi-argument forward rule
function ruleSPNonlinearEOutNGX(g::Function, # Needs to be in front of Vararg
                                msg_out::Nothing,
                                msgs_in::Vararg{Message{<:Gaussian, V}}) where V<:VariateType

    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearization(V, g, ms_fw_in)
    (m_fw_in, V_fw_in, _) = concatenateGaussianMV(ms_fw_in, Vs_fw_in)

    return Message(V, GaussianMeanVariance, m=A*m_fw_in + b, v=A*V_fw_in*A')
end

# Backward rule with given inverse
function ruleSPNonlinearEIn1GG(g::Function,
                               g_inv::Function,
                               msg_out::Message{F, V},
                               msg_in1::Nothing) where {F<:Gaussian, V<:VariateType}

    (m_out, V_out) = unsafeMeanCov(msg_out.dist)
    (A, b) = localLinearization(V, g_inv, m_out)

    return Message(V, GaussianMeanVariance, m=A*m_out + b, v=A*V_out*A')
end

# Multi-argument backward rule with given inverse
function ruleSPNonlinearEInGX(g::Function, # Needs to be in front of Vararg
                              g_inv::Function,
                              msg_out::Message{<:Gaussian, V},
                              msgs_in::Vararg{Union{Message{<:Gaussian, V}, Nothing}}) where V<:VariateType

    (ms, Vs) = collectStatistics(msg_out, msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearization(V, g_inv, ms)
    (mc, Vc) = concatenateGaussianMV(ms, Vs)

    return Message(V, GaussianMeanVariance, m=A*mc, v=A*Vc*A')
end

# Backward rule with unknown inverse
function ruleSPNonlinearEIn1GG(g::Function,
                               msg_out::Message{F1, V},
                               msg_in1::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:VariateType}

    m_in1 = unsafeMean(msg_in1.dist)
    d_out = convert(ProbabilityDistribution{V, GaussianMeanPrecision}, msg_out.dist)
    m_out = d_out.params[:m]
    W_out = d_out.params[:w]
    (A, b) = localLinearization(V, g, m_in1)

    return Message(V, GaussianWeightedMeanPrecision, xi=A'*W_out*(m_out - b), w=A'*W_out*A)
end

# Multi-argument backward rule with unknown inverse
function ruleSPNonlinearEInGX(g::Function,
                              inx::Int64, # Index of inbound interface inx
                              msg_out::Message{<:Gaussian, V},
                              msgs_in::Vararg{Message{<:Gaussian, V}}) where V<:VariateType

    # Approximate joint inbounds
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearization(V, g, ms_fw_in)
    
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in)
    m_fw_out = A*m_fw_in + b
    V_fw_out = A*V_fw_in*A'
    C_fw = V_fw_in*A'

    # RTS Smoother
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    (m_in, V_in) = smoothRTS(m_fw_out, V_fw_out, C_fw, m_fw_in, V_fw_in, m_bw_out, V_bw_out)

    # Marginalize joint belief on in's
    (m_inx, V_inx) = marginalizeGaussianMV(V, m_in, V_in, ds, inx) # Marginalization is overloaded on VariateType V
    W_inx = cholinv(V_inx) # Convert to canonical statistics
    xi_inx = W_inx*m_inx

    # Divide marginal on inx by forward message
    (xi_fw_inx, W_fw_inx) = unsafeWeightedMeanPrecision(msgs_in[inx].dist)
    xi_bw_inx = xi_inx - xi_fw_inx
    W_bw_inx = W_inx - W_fw_inx # Note: subtraction might lead to posdef inconsistencies

    return Message(V, GaussianWeightedMeanPrecision, xi=xi_bw_inx, w=W_bw_inx)
end

function ruleMNonlinearEInGX(g::Function,
                             msg_out::Message{<:Gaussian, V},
                             msgs_in::Vararg{Message{<:Gaussian, V}}) where V<:VariateType

    # Approximate joint inbounds
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearization(V, g, ms_fw_in)
    
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in)
    m_fw_out = A*m_fw_in + b
    V_fw_out = A*V_fw_in*A'
    C_fw = V_fw_in*A'

    # RTS Smoother
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    (m_in, V_in) = smoothRTS(m_fw_out, V_fw_out, C_fw, m_fw_in, V_fw_in, m_bw_out, V_bw_out)

    return ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=m_in, v=V_in)
end
