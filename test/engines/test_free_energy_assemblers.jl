module FreeEnergyAssemblersTest

using Test
using ForneyLab
using ForneyLab: assembleCountingNumbers!, setTargets!, assembleFreeEnergy!

@testset "assembleCountingNumbers!" begin
    fg = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(0.0, 1.0)
    @RV z = x + y
    @RV w ~ GaussianMeanPrecision(z, 1.0)
    placeholder(w, :w)

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, free_energy=true)
    assembleCountingNumbers!(pfz)

    cl = first(pf.target_clusters)
    nd1 = fg.nodes[:gaussianmeanprecision_1]
    nd2 = fg.nodes[:gaussianmeanprecision_2]
    nd3 = fg.nodes[:gaussianmeanprecision_3]

    @test pfz.entropy_counting_numbers[x] == 0
    @test pfz.entropy_counting_numbers[cl] == 1
    @test pfz.energy_counting_numbers[nd1] == 1
    @test pfz.energy_counting_numbers[nd2] == 1
    @test pfz.energy_counting_numbers[nd3] == 1

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, free_energy=false)
    @test_throws Exception assembleCountingNumbers!(pfz)

    # Regression test
    fg = FactorGraph()
    @RV x ~ GaussianMeanVariance(placeholder(:mx), placeholder(:vx))
    @RV y ~ GaussianMeanVariance(placeholder(:my), placeholder(:vy))
    @RV z = dot(x, 1.0) + dot(y, 1.0)
    placeholder(z, :z)

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    setTargets!(pf, pfz, free_energy=true)
    assembleCountingNumbers!(pfz)

    cl = first(pf.target_clusters)
    nd1 = fg.nodes[:gaussianmeanvariance_1]
    nd2 = fg.nodes[:gaussianmeanvariance_2]
    v1 = fg.variables[:variable_1]
    v2 = fg.variables[:variable_2]

    @test pfz.entropy_counting_numbers[cl] == 1
    @test pfz.entropy_counting_numbers[v1] == -1
    @test pfz.entropy_counting_numbers[v2] == -1
    @test pfz.entropy_counting_numbers[x] == 1
    @test pfz.entropy_counting_numbers[y] == 1
    @test pfz.energy_counting_numbers[nd1] == 1
    @test pfz.energy_counting_numbers[nd2] == 1
end

@testset "assembleFreeEnergy!" begin
    fg = FactorGraph()
    @RV x ~ GaussianMeanPrecision(0.0, 1.0)
    @RV y ~ GaussianMeanPrecision(0.0, 1.0)
    @RV z = x + y
    @RV w ~ GaussianMeanPrecision(z, 1.0)
    placeholder(w, :w)

    pfz = PosteriorFactorization()
    pf = PosteriorFactor(fg)
    algo = messagePassingAlgorithm(Variable[], free_energy=true)
    assembleFreeEnergy!(algo)

    cl = first(pf.target_clusters)
    nd1 = fg.nodes[:gaussianmeanprecision_1]
    nd2 = fg.nodes[:gaussianmeanprecision_2]
    nd3 = fg.nodes[:gaussianmeanprecision_3]

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
