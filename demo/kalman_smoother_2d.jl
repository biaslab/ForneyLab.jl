
# RTS smoothing for a linear Gaussian state space model.
# based on https://github.com/probml/pmtk3/blob/master/demos/kalmanTrackingDemo.m
# We consider the example from "Machine learning: a probabilistic perspective"
# (Murphy, 2012) fig 18.1.

using ForneyLab,  LinearAlgebra, Test
import Distributions # "Multivariate" clashes with ForneyLab
#using PDMats

function make_diagonal(v)
    d = length(v)
    A = zeros(d,d)
    A[diagind(A)] = v
    return A
end

function make_params()
    F = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1.0];
    H = [1.0 0 0 0; 0 1 0 0];
    nobs, nhidden = size(H)
    #Q = PDiagMat(fill(0.001, nhidden))
    #Q = make_diagonal(fill(0.001, nhidden))
    e = 0.001; Q = [e 0 0 0; 0 e 0 0; 0 0 e 0; 0 0 0 e]
    #R = PDiagMat(fill(1.0, nobs))
    #R = make_diagonal(fill(1.0, nobs))
    e = 1.0; R = [e 0; 0 e]
    mu0 = [8, 10, 1, 0.0];
    #V0 = PDiagMat(fill(1, nhidden))
    #V0 = make_diagonal(fill(1, nhidden))
    e = 1.0; V0 = [e 0 0 0; 0 e 0 0; 0 0 e 0; 0 0 0 e]
    params = (mu0 = mu0, V0 = V0, F = F, H = H, Q = Q, R = R)
    return params
end

function lin_gauss_ssm_sample(T, params)
    # z(t+1) = F * z(t-1) + N(0,Q) so F is D*D
    # y(t) = H * z(t) + N(0, R) so H is O*D
    # Returns zs: H*T, ys: O*T
    F = params.F; H = params.H; Q = params.Q; R = params.R; mu0 = params.mu0; V0 = params.V0;
    nobs, nhidden = size(H)
    zs = Array{Float64}(undef, nhidden, T)
    ys = Array{Float64}(undef, nobs, T)
    prior_z = Distributions.MvNormal(mu0, V0)
    process_noise_dist = Distributions.MvNormal(Q)
    obs_noise_dist = Distributions.MvNormal(R)
    zs[:,1] = rand(prior_z)
    ys[:,1] = H*zs[:,1] + rand(obs_noise_dist)
    for t=2:T
        zs[:,t] = F*zs[:,t-1] + rand(process_noise_dist)
        ys[:,t] = H*zs[:,t] + rand(obs_noise_dist)
    end
    return zs, ys
end

function convert_2d_matrix_to_vec_of_vec(M::Array{Float64,2})
    N, T = size(M)
    v = Vector{Vector{Float64}}(undef, T)
    for i=1:T
        v[i] = M[:,i]
    end
    return v
end

function make_data(T, params)
    nobs, nhidden = size(params.H)
    (zs, ys) = lin_gauss_ssm_sample(T, params)
    #zs = randn(nhidden, T)
    #ys = randn(nobs, T)
    # ForneyLab cannot handle matrix observations
    # https://github.com/biaslab/ForneyLab.jl/issues/17
    y_data =  convert_2d_matrix_to_vec_of_vec(ys)
end


function make_factor_graph(T, params)
    F = params.F; H = params.H; Q = params.Q; R = params.R; mu0 = params.mu0; V0 = params.V0;
    g = FactorGraph()
    @RV x0 ~ GaussianMeanVariance(mu0, V0)
    x = Vector{Variable}(undef, T)
    y = Vector{Variable}(undef, T)
    x_t_prev = x0
    for t = 1:T
        #global x_t_prev
        mu_x = F * x_t_prev
        @RV [id=:x_*t] x[t] ~ GaussianMeanVariance(mu_x, Q)
        #@RV x[t] ~ GaussianMeanVariance(mu_x, Q)
        #x[t].id = Symbol("x_", t) #:x_t;
        mu_y = H * x[t]
        @RV y[t] ~ GaussianMeanVariance(mu_y, R)
        nobs, nhidden = size(H)
        placeholder(y[t], :y, dims=(nobs,), index=t)
        x_t_prev = x[t]
    end
    fg = (x0 = x0, x = x, y = y)
    return fg
end


function make_fg_inference(fg)
    println("generating code")
    algo = Meta.parse(sumProductAlgorithm(fg.x))
    eval(algo) # Compile the step! function
    function infer(y_data)
        T = length(y_data)
        D = length(y_data[1])
        id = Symbol("y")
        data = Dict(id => y_data)
        #data = Dict(fg.y => y_data)
        println("running forwards backwards")
        marginals = step!(data)
        println("computing marginals")
        m_x = [mean(marginals[:x_*t]) for t = 1:T] # T of D arrays
        V_x = [cov(marginals[:x_*t]) for t = 1:T] # T of DxD arrays
        return m_x, V_x
    end
    return infer
end

 # Source:
# https://github.com/QuantEcon/QuantEcon.lectures.code/blob/master/kalman/gaussian_contours.jl
function bivariate_normal(X::Matrix, Y::Matrix, σ_x::Real=1.0, σ_y::Real=1.0,
                          μ_x::Real=0.0, μ_y::Real=0.0, σ_xy::Real=0.0)
    Xμ = X .- μ_x
    Yμ = Y .- μ_y
    ρ = σ_xy/(σ_x*σ_y)
    z = Xμ.^2/σ_x^2 + Yμ.^2/σ_y^2 - 2*ρ.*Xμ.*Yμ/(σ_x*σ_y)
    denom = 2π*σ_x*σ_y*sqrt(1-ρ^2)
    return exp.(-z/(2*(1-ρ^2))) ./ denom
end


# https://github.com/probml/pmtk3/blob/master/matlabTools/graphics/gaussPlot2d.m
function plot_gauss2d(m, C)
    U = eigvecs(C)
    D = eigvals(C)
    N = 100
    t = range(0, stop=2*pi, length=N)
    xy = zeros(Float64, 2, N)
    xy[1,:] = cos.(t)
    xy[2,:] = sin.(t)
    k = sqrt(6) # approx sqrt(chi2inv(0.95, 2)) = 2.45
    k = 0.5 # random hack to make ellipses smaller
    w = (k * U * Diagonal(sqrt.(D))) * xy # 2*N
    Plots.scatter!([m[1]], [m[2]], marker=:star, label="")
    handle = plot!(w[1,:] .+ m[1], w[2,:] .+ m[2], label="")
    return handle
end

function plot_gauss2d_test()
    m = [0.0, 0.0]
    #C = randn(2,2); C = C*C';
    C = [1.0 0.0; 0.0 3.0];
    @test isposdef(C)
    Plots.scatter([m[1]], [m[2]], marker=:star)
    plot_2dgauss(m, C)
end


function compute_gauss_marginals(m_x, V_x, keep)
    mm_x = [m_x[t][keep] for t=1:T]
    VV_x = [V_x[t][keep, keep] for t=1:T]
    return mm_x, VV_x
end


Random.seed!(1)
T = 10
model = make_params() # cannot use "params" because of "Flux.params"
y_data = make_data(T, model) # T vector of D vectors
ys = hcat(y_data...) # D*T matrix
fg = make_factor_graph(T, model)
infer = make_fg_inference(fg)
mz, Vz = infer(y_data)
mz2, Vz2 = compute_gauss_marginals(mz, Vz, 1:2)

using Plots; pyplot()
closeall()
plt = scatter(ys[1,:], ys[2,:], label="obs", reuse=false)
xlims!(minimum(ys[1,:])-1, maximum(ys[1,:])+1)
ylims!(minimum(ys[2,:])-1, maximum(ys[2,:])+1)
display(plt)
for t=1:T
    plt = plot_gauss2d(mz2[t], Vz2[t])
    #plt = annotate!(mz2[t][1], mz2[t][2], "t=$t");
    display(plt)
end
