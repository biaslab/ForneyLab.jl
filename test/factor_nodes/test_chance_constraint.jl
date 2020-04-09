module ChanceConstraintTest

using Test
using ForneyLab
using ForneyLab: outboundType, isApplicable
using ForneyLab: SPChanceConstraintOutG
u
g(x::Float64) = 1.0*(x > 0.0)


#-------------
# Update rules
#-------------

@testset "SPChanceConstraintOutG" begin
    @test true == false
end

end # module