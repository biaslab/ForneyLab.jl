using ForneyLab

GaussianMixtureNodePar(id=:gm)
PriorNode(Categorical([0.25,0.25,0.5]),id=:z)
PriorNode(Dirichlet([10.,10.,20.]),id=:pi)
PriorNode(PartitionedDistribution([Gaussian(m=5., V=50.),Gaussian(m=10., V=50.),Gaussian(m=20., V=100.)]),id=:m)
PriorNode(PartitionedDistribution([Gamma(10.,10.),Gamma(2.,2.),Gamma(2.,5.)]),id=:w)
TerminalNode(Mixture([Gaussian(m=0.0, V=huge), Gaussian(m=0.0, V=huge), Gaussian(m=0.0, V=huge)],[0.5,0.5]),id=:x)

Edge(n(:z).i[:out],n(:gm).i[:z],id=:e_z)
Edge(n(:m).i[:out],n(:gm).i[:m], id=:e_m)
Edge(n(:pi).i[:out],n(:gm).i[:pi],id=:e_pi)
Edge(n(:w).i[:out],n(:gm).i[:w],id=:e_w)
Edge(n(:x).i[:out], n(:gm).i[:x],id=:e_x)

x_est=attachWriteBuffer(n(:x).i[:out].partner)

algo=SumProduct()
run(algo)

println(x_est)
