@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Beta},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPBetaOutNMM)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Nothing, Message),
                :name          => SPBetaMNM)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Message, Nothing),
                :name          => SPBetaMMN)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBBetaOut)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBBetaA)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBBetaB)


mutable struct SPBetaOutMCNMM <: SumProductRule{Beta} end
outboundType(::Type{SPBetaOutMCNMM}) = Message{SampleList}
function isApplicable(::Type{SPBetaOutMCNMM}, input_types::Vector{Type})
  nothing_inputs = 0
  point_inputs = 0

  for input_type in input_types
      if input_type == Nothing
          nothing_inputs += 1
      elseif matches(input_type, Message{PointMass})
          point_inputs += 1
      end
  end

  return (nothing_inputs == 1) && (point_inputs != 2)
end
