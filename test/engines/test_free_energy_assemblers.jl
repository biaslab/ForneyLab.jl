module FreeEnergyAssemblersTest

using Test
using ForneyLab
using ForneyLab: assembleCountingNumbers!, setTargets!, assembleFreeEnergy!

@testset "assembleCountingNumbers!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(0.0, 1.0)
    @RV z = x + y
    @RV w ~ GaussianMeanPrecision(z, 1.0)
    placeholder(w, :w)

    algo = Algorithm()
    rf = RecognitionFactor(algo)
    setTargets!(rf, algo, free_energy=true)
    assembleCountingNumbers!(algo)

    cl = first(rf.clusters)
    nd1 = g.nodes[:gaussianmeanprecision_1]
    nd2 = g.nodes[:gaussianmeanprecision_2]
    nd3 = g.nodes[:gaussianmeanprecision_3]

    @test algo.entropy_counting_numbers[x] == 0
    @test algo.entropy_counting_numbers[cl] == 1
    @test algo.energy_counting_numbers[nd1] == 1
    @test algo.energy_counting_numbers[nd2] == 1
    @test algo.energy_counting_numbers[nd3] == 1
end

@testset "assembleFreeEnergy!" begin
    g = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(0.0, 1.0)
    @RV z = x + y
    @RV w ~ GaussianMeanPrecision(z, 1.0)
    placeholder(w, :w)

    algo = Algorithm()
    sumProductAlgorithm(Variable[], algo, free_energy=true)
    rf = algo.recognition_factors[Symbol("")]

    cl = first(rf.clusters)
    nd1 = g.nodes[:gaussianmeanprecision_1]
    nd2 = g.nodes[:gaussianmeanprecision_2]
    nd3 = g.nodes[:gaussianmeanprecision_3]

    @test algo.entropies[1][:counting_number] == 1
    @test algo.entropies[1][:inbound].target == cl
    @test algo.average_energies[1][:counting_number] == 1
    @test length(algo.average_energies[1][:inbounds]) == 3
    @test algo.average_energies[2][:counting_number] == 1
    @test length(algo.average_energies[2][:inbounds]) == 3
    @test algo.average_energies[3][:counting_number] == 1
    @test length(algo.average_energies[3][:inbounds]) == 3
end

end # module
