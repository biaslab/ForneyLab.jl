# Bayes rule for a linear Gaussian model
# X -> Y, observe Y, infer X
# Y ~ N(A x + b, Sigma_y)
# X ~ N(mu_x, Sigma_x)

# We consider the example from "Machine learning: a probabilistic perspective"
# (Murphy, 2012) sec 4.4.2.2.
# where A=I, b=0, mu=[0.5 0.5], Sigma_x=0.1I, Sigma_y=0.1*[2 1; 1 1]
# Corresponding Matlab code is here:
# https://github.com/probml/pmtk3/blob/master/demos/gaussInferParamsMean2d.m

using ForneyLab
using LinearAlgebra, Distributions, Test

function make_data(params, n_data)
    y_mean = params.A * params.true_x  + params.b
    y_dist = Distributions.MvNormal(y_mean, params.Sigma_y)
    y_data_all = rand(y_dist, n_data) # n_data x obs_dim
    return y_data_all
end

function make_factor_graph(params)
    g = FactorGraph()
    @RV x ~ GaussianMeanVariance(params.mu_x, params.Sigma_x)
    mu_y = params.A * x + params.b
    @RV y ~ GaussianMeanVariance(mu_y, params.Sigma_y)
    placeholder(y, :y, dims=(length(params.b),))
    fg = (x = x, y = y)
    return fg
end

function make_fg_inference(fg)
    algo = Meta.parse(sumProductAlgorithm(fg.x))
    eval(algo) # compile the step! function
    function infer(y_data)
        data = Dict(fg.y.id => y_data)
        marginals = step!(data);
        post = marginals[:x]
        post_gauss = Distributions.MvNormal(mean(post), cov(post))
        return post_gauss
    end
    return infer
end

function bayes_rule_lin_gauss(params, y)
    # Implements eqn 4.125 of MLAPP
    A = params.A
    At = permutedims(A)
    b = params.b
    prior_prec = inv(params.Sigma_x)
    prior_mean = params.mu_x
    obs_prec = inv(params.Sigma_y)
    post_prec = prior_prec + At * obs_prec * A
    post_cov = inv(post_prec)
    post_mean = post_cov*(At * obs_prec * (y-b) + prior_prec * prior_mean)
    post_gauss = Distributions.MvNormal(post_mean, post_cov)
    return post_gauss
end

#=
params = Dict(
    :true_x => [0.5, 0.5],
     :mu_x => [0.0, 0.0],
     :Sigma_x => 0.1 * eye(2),
     :A => eye(2),
     :b => [0.0, 0.0],
     :Sigma_y => 0.1 .* [2 1; 1 1]
     )
=#

 params_gen =
     (true_x = [0.5, 0.5],
      mu_x = [0.0, 0.0],
      Sigma_x = 0.1 * eye(2),
      A = eye(2),
      b = [0.0, 0.0],
      Sigma_y = 0.1 .* [2 1; 1 1]
      )

n_data = 10
Random.seed!(1)
y_data_all = make_data(params_gen, n_data)
y_data_mean = vec(mean(y_data_all, dims=2))

params =
    (true_x = [0.5, 0.5],
     mu_x = [0.0, 0.0],
     Sigma_x = 0.1 * eye(2),
     A = eye(2),
     b = [0.0, 0.0],
     Sigma_y = 0.1 .* [2 1; 1 1] / n_data # since using average of observations
     )
fg = make_factor_graph(params)
infer = make_fg_inference(fg)
post_gauss = infer(y_data_mean)
post_gauss2 = bayes_rule_lin_gauss(params, y_data_mean)
@test isapprox(mean(post_gauss), mean(post_gauss2))
@test isapprox(cov(post_gauss), cov(post_gauss2))

using Plots; pyplot()
closeall()
xrange = -1.5:0.01:1.5
yrange = xrange
plt = scatter(y_data_all[1,:], y_data_all[2,:], label="obs", reuse=false)
scatter!([params.true_x[1]], [params.true_x[2]], label="truth")
xlims!(minimum(xrange), maximum(xrange))
ylims!(minimum(xrange), maximum(xrange))
title!("Data")
savefig("Figures/gauss-bayes-rule-2d-data.png")
display(plt)

prior_gauss = Distributions.MvNormal(params.mu_x, params.Sigma_x)
plt = contour(xrange, yrange, (x,y)->Distributions.pdf(prior_gauss,[x,y]),
    reuse=false, title="prior")
savefig("Figures/gauss-bayes-rule-2d-prior.png")
display(plt)

plt = contour(xrange, yrange, (x,y)->Distributions.pdf(post_gauss,[x,y]),
    reuse=false, title = "Posterior")
savefig("Figures/gauss-bayes-rule-2d-post.png")
display(plt)

# foo
