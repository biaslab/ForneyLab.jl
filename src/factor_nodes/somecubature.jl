const product = Iterators.product
const PIterator = Iterators.ProductIterator

function generate_multidim_points(n::Int, p::Int)
   sigma_points, sigma_weights = gausshermite(p)
   points_iter = product(repeat([sigma_points],n)...)
   weights_iter = product(repeat([sigma_weights],n)...)
   return points_iter, weights_iter
end
function gaussHermiteCubature(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::PIterator, weights_iter::PIterator)
   result = 0.0
   # println(d.params[:v])
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result += weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
function gaussHermiteCubature1D(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::Array, weights_iter::Array)
   result = 0.0
   # println(d.params[:v])
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result += weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
function gaussHermiteQuadrature(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::Array, weights_iter::Array)
   result = 0.0
   # println(d.params[:v])
   std = sqrt(d.params[:v])
   for (point, weight) in zip(points_iter, weights_iter)
       result += weight*g(d.params[:m] + std*point)
   end
   return result
end
function gaussHermiteCubatureMean(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::PIterator, weights_iter::PIterator)
   result = zeros(dims(d))
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result = result + weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
function gaussHermiteCubatureCov(g::Function, d::ProbabilityDistribution{Univariate,GaussianMeanVariance}, points_iter::PIterator, weights_iter::PIterator)
   result = zeros(dims(d),dims(d))
   sqrtP = sqrt(d.params[:v])
   for (point_tuple, weights) in zip(points_iter, weights_iter)
       weight = prod(weights)
       point = collect(point_tuple)
       result = result + weight.*g(d.params[:m] + sqrtP*point)
   end
   return result
end
