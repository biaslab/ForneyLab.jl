module CompositeTest

using Base.Test
using ForneyLab
import ForneyLab: @composite

@testset "@composite" begin
    @test true == false
end


#-------------
# Update rules
#-------------

@testset "Custom SumProductRule for Composite" begin
    @test true == false
end

end #module