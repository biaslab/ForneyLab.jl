export ruleVBWishartOut

ruleVBWishartOut(   dist_out::Any,
                    dist_v::ProbabilityDistribution{MatrixVariate},
                    dist_nu::ProbabilityDistribution{Univariate}) =
    Message(MatrixVariate(Wishart, v=unsafeMean(dist_v), nu=unsafeMean(dist_nu)))