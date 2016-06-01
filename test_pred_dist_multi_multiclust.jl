using ForneyLab

GaussianMixtureNode(id=:gm)
PriorNode(Categorical([0.25,0.25,0.5]),id=:z)
PriorNode(Dirichlet([10.,10.,20.]),id=:pi)
PriorNode(PartitionedDistribution([MvGaussian(m=[5., 10.], V=[50. 0.; 0. 50.]),MvGaussian(m=[10., 5.], V=[50. 0.; 0. 50.]),MvGaussian(m=[20., 10.], V=[1. 0.; 0. 1.])]),id=:m)
PriorNode(PartitionedDistribution([Wishart(V=[10. 0.;0. 10.],nu=10.),Wishart(V=[10. 0.;0. 10.],nu=10.),Wishart(V=[10. 0.;0. 10.],nu=10.)]),id=:w)
TerminalNode(Mixture([MvGaussian(m=[0.0; 0.0], V=[huge 0; 0 huge]), MvGaussian(m=[0.0; 0.0], V=[huge 0; 0 huge]), MvGaussian(m=[0.0; 0.0], V=[huge 0; 0 huge])],[1./3.,1./3.,1./3.]),id=:x)

Edge(n(:z).i[:out],n(:gm).i[:z],id=:e_z)
Edge(n(:m).i[:out],n(:gm).i[:m], id=:e_m)
Edge(n(:pi).i[:out],n(:gm).i[:pi],id=:e_pi)
Edge(n(:w).i[:out],n(:gm).i[:w],id=:e_w)
Edge(n(:x).i[:out], n(:gm).i[:x],id=:e_x)

x_est=attachWriteBuffer(n(:x).i[:out].partner)

algo=SumProduct()
run(algo)

println(x_est)
