@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Beta},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPBetaOutNPP)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Nothing, Message),
                :name          => SPBetaAMNM)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Message, Nothing),
                :name          => SPBetaBMMN)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBBetaOut)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Function},
                      :inbound_types => (Distribution, Nothing, Distribution),
                      :name          => VBBetaA)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Function},
                      :inbound_types => (Distribution, Distribution, Nothing),
                      :name          => VBBetaB)


mutable struct SPBetaOutNMM <: SumProductRule{Beta} end
outboundType(::Type{SPBetaOutNMM}) = Message{SampleList}
function isApplicable(::Type{SPBetaOutNMM}, input_types::Vector{Type})
  nothing_inputs = 0
  point_inputs = 0

  for input_type in input_types
      if input_type == Nothing
          nothing_inputs += 1
      elseif input_type << Message{PointMass}
          point_inputs += 1
      end
  end

  return (nothing_inputs == 1) && (point_inputs != 2)
end
