import Base: vec

export
ruleSPDeltaEOutNG,
ruleSPDeltaEOutNGX,
ruleSPDeltaEIn1GG,
ruleSPDeltaEInGX,
ruleMDeltaEInGX

"""
Concatenate a vector (of vectors and floats) and return with original dimensions (for splitting)
"""
function concatenate(xs::Vector)
    ds = size.(xs) # Extract dimensions
    x = vcat(vec.(xs)...)

    return (x, ds)
end

ForneyLab.vec(d::Float64) = [d] # Extend vectorization to Float

"""
Return local linearization of g around expansion point x_hat
for Delta node with single input interface
"""
function localLinearizationSingleIn(g::Function, x_hat::Float64)
    a = ForwardDiff.derivative(g, x_hat)
    b = g(x_hat) - a*x_hat

    return (a, b)
end

function localLinearizationSingleIn(g::Function, x_hat::Vector{Float64})
    A = ForwardDiff.jacobian(g, x_hat)
    b = g(x_hat) - A*x_hat

    return (A, b)
end

"""
Return local linearization of g around expansion point x_hat
for Delta node with multiple input interfaces
"""
function localLinearizationMultiIn(g::Function, x_hat::Vector{Float64})
    g_unpacked(x::Vector) = g(x...)
    A = ForwardDiff.gradient(g_unpacked, x_hat)'
    b = g(x_hat...) - A*x_hat

    return (A, b)
end

function localLinearizationMultiIn(g::Function, x_hat::Vector{Vector{Float64}})
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
function ruleSPDeltaEOutNG(g::Function,
                               msg_out::Nothing,
                               msg_in1::Message{<:Gaussian})
    
    (m_in1, V_in1) = unsafeMeanCov(msg_in1.dist)
    (A, b) = localLinearizationSingleIn(g, m_in1)
    m = A*m_in1 + b
    V = A*V_in1*A'

    return Message(variateType(m), GaussianMeanVariance, m=m, v=V)
end

# Multi-argument forward rule
function ruleSPDeltaEOutNGX(g::Function, # Needs to be in front of Vararg
                                msg_out::Nothing,
                                msgs_in::Vararg{Message{<:Gaussian}})

    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearizationMultiIn(g, ms_fw_in)
    (m_fw_in, V_fw_in, _) = concatenateGaussianMV(ms_fw_in, Vs_fw_in)
    m = A*m_fw_in + b
    V = A*V_fw_in*A'

    return Message(variateType(m), GaussianMeanVariance, m=m, v=V)
end

# Backward rule with given inverse
function ruleSPDeltaEIn1GG(g::Function,
                               g_inv::Function,
                               msg_out::Message{<:Gaussian},
                               msg_in1::Nothing)

    (m_out, V_out) = unsafeMeanCov(msg_out.dist)
    (A, b) = localLinearizationSingleIn(g_inv, m_out)
    m = A*m_out + b
    V = A*V_out*A'

    return Message(variateType(m), GaussianMeanVariance, m=m, v=V)
end

# Multi-argument backward rule with given inverse
function ruleSPDeltaEInGX(g::Function, # Needs to be in front of Vararg
                              g_inv::Function,
                              msg_out::Message{<:Gaussian},
                              msgs_in::Vararg{Union{Message{<:Gaussian}, Nothing}})

    (ms, Vs) = collectStatistics(msg_out, msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearizationMultiIn(g_inv, ms)
    (mc, Vc) = concatenateGaussianMV(ms, Vs)
    m = A*mc + b
    V = A*Vc*A'

    return Message(variateType(m), GaussianMeanVariance, m=m, v=V)
end

# Backward rule with unknown inverse
function ruleSPDeltaEIn1GG(g::Function,
                               msg_out::Message{<:Gaussian},
                               msg_in1::Message{<:Gaussian})

    m_in1 = unsafeMean(msg_in1.dist)
    (m_out, W_out) = unsafeMeanPrecision(msg_out.dist)
    (A, b) = localLinearizationSingleIn(g, m_in1)
    xi = A'*W_out*(m_out - b)
    W = A'*W_out*A

    return Message(variateType(xi), GaussianWeightedMeanPrecision, xi=xi, w=W)
end

# Multi-argument backward rule with unknown inverse
function ruleSPDeltaEInGX(g::Function,
                              inx::Int64, # Index of inbound interface inx
                              msg_out::Message{<:Gaussian},
                              msgs_in::Vararg{Message{<:Gaussian}})

    # Approximate joint inbounds
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearizationMultiIn(g, ms_fw_in)
    
    (m_fw_in, V_fw_in, ds) = concatenateGaussianMV(ms_fw_in, Vs_fw_in)
    m_fw_out = A*m_fw_in + b
    V_fw_out = A*V_fw_in*A'
    C_fw = V_fw_in*A'

    # RTS Smoother
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    (m_in, V_in) = smoothRTS(m_fw_out, V_fw_out, C_fw, m_fw_in, V_fw_in, m_bw_out, V_bw_out)

    # Marginalize joint belief on in's
    (m_inx, V_inx) = marginalizeGaussianMV(m_in, V_in, ds, inx)
    W_inx = cholinv(V_inx) # Convert to canonical statistics
    xi_inx = W_inx*m_inx

    # Divide marginal on inx by forward message
    (xi_fw_inx, W_fw_inx) = unsafeWeightedMeanPrecision(msgs_in[inx].dist)
    xi_bw_inx = xi_inx - xi_fw_inx
    W_bw_inx = W_inx - W_fw_inx # Note: subtraction might lead to posdef violations

    return Message(variateType(xi_bw_inx), GaussianWeightedMeanPrecision, xi=xi_bw_inx, w=W_bw_inx)
end

function ruleMDeltaEInGX(g::Function,
                             msg_out::Message{<:Gaussian},
                             msgs_in::Vararg{Message{<:Gaussian}})

    # Approximate joint inbounds
    (ms_fw_in, Vs_fw_in) = collectStatistics(msgs_in...) # Returns arrays with individual means and covariances
    (A, b) = localLinearizationMultiIn(g, ms_fw_in)
    
    (m_fw_in, V_fw_in, _) = concatenateGaussianMV(ms_fw_in, Vs_fw_in)
    m_fw_out = A*m_fw_in + b
    V_fw_out = A*V_fw_in*A'
    C_fw = V_fw_in*A'

    # RTS Smoother
    (m_bw_out, V_bw_out) = unsafeMeanCov(msg_out.dist)
    (m_in, V_in) = smoothRTS(m_fw_out, V_fw_out, C_fw, m_fw_in, V_fw_in, m_bw_out, V_bw_out)

    return Distribution(Multivariate, GaussianMeanVariance, m=m_in, v=V_in)
end
