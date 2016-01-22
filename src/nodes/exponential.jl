############################################
# ExponentialNode
############################################
# Description:
#   Maps a Gaussian to a log-normal distribution.
#   Derivations can be found in the derivations document.
#
#    in         out
#   ----->[exp]----->
#
#   f(in,out) = δ(out - exp(in))
#
# Interfaces:
#   1 i[:in], 2 i[:out]
#
# Construction:
#   ExponentialNode(id=:my_node)
#
############################################

export ExponentialNode

type ExponentialNode <: Node
    id::Symbol
    interfaces::Array{Interface,1}
    i::Dict{Symbol,Interface}

    function ExponentialNode(; id=generateNodeId(ExponentialNode))
        self = new(id, Array(Interface, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        for (iface_index, iface_handle) in enumerate([:in, :out])
            self.i[iface_handle] = self.interfaces[iface_index] = Interface(self)
        end

        return self
    end
end

isDeterministic(::ExponentialNode) = true


############################################
# Gaussian update functions
############################################

function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Type{Val{2}},
                        outbound_dist::LogNormalDistribution,
                        msg_in::Message{GaussianDistribution},
                        msg_out::Any)

    # TODO: unrepress
    # isProper(msg_in.payload) || error("Improper input distributions are not supported")
    ensureParameters!(msg_in.payload, (:m, :V))

    outbound_dist.m = msg_in.payload.m
    outbound_dist.s = msg_in.payload.V

    return outbound_dist
end

function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Type{Val{1}},
                        outbound_dist::GaussianDistribution,
                        msg_in::Any,
                        msg_out::Message{LogNormalDistribution})

    # TODO: unrepress
    # isProper(msg_out.payload) || error("Improper input distributions are not supported")

    outbound_dist.m = msg_out.payload.m
    outbound_dist.V = msg_out.payload.s
    outbound_dist.W = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end


############################################
# DeltaDistribution update functions
############################################

function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Type{Val{2}},
                        outbound_dist::DeltaDistribution{Float64},
                        msg_in::Message{DeltaDistribution{Float64}},
                        msg_out::Any)

    outbound_dist.m = exp(msg_in.payload.m)
    return outbound_dist
end

function sumProduct!(   node::ExponentialNode,
                        outbound_interface_index::Type{Val{1}},
                        outbound_dist::DeltaDistribution{Float64},
                        msg_in::Any,
                        msg_out::Message{DeltaDistribution{Float64}})

    outbound_dist.m = log(msg_out.payload.m)
    return outbound_dist
end


############################################
# MvGaussian update functions
############################################

function sumProduct!{dims<:Int64}(  node::ExponentialNode,
                                    outbound_interface_index::Type{Val{2}},
                                    outbound_dist::MvLogNormalDistribution{dims},
                                    msg_in::Message{MvGaussianDistribution{dims}},
                                    msg_out::Any)

    # TODO: unrepress
    # isProper(msg_in.payload) || error("Improper input distributions are not supported")
    ensureParameters!(msg_in.payload, (:m, :V))

    outbound_dist.m = deepcopy(msg_in.payload.m)
    outbound_dist.S = deepcopy(msg_in.payload.V)

    return outbound_dist
end

function sumProduct!{dims<:Int64}(  node::ExponentialNode,
                                    outbound_interface_index::Type{Val{1}},
                                    outbound_dist::MvGaussianDistribution{dims},
                                    msg_in::Any,
                                    msg_out::Message{MvLogNormalDistribution{dims}})

    # TODO: unrepress
    # isProper(msg_out.payload) || error("Improper input distributions are not supported")

    outbound_dist.m = deepcopy(msg_out.payload.m)
    outbound_dist.V = deepcopy(msg_out.payload.S)
    outbound_dist.W = NaN
    outbound_dist.xi = NaN

    return outbound_dist
end


############################################
# MvDeltaDistribution update functions
############################################

function sumProduct!{T<:MvDeltaDistribution{Float64}}(  node::ExponentialNode,
                                                        outbound_interface_index::Type{Val{2}},
                                                        outbound_dist::T,
                                                        msg_in::Message{T},
                                                        msg_out::Any)

    outbound_dist.m = exp(msg_in.payload.m)
    return outbound_dist
end

function sumProduct!{T<:MvDeltaDistribution{Float64}}(  node::ExponentialNode,
                                                        outbound_interface_index::Type{Val{1}},
                                                        outbound_dist::T,
                                                        msg_in::Any,
                                                        msg_out::Message{T})

    outbound_dist.m = log(msg_out.payload.m)
    return outbound_dist
end
