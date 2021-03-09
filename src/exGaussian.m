function f = exGaussian( x, mu, sigma, lambda )
% EXGAUSSIAN implements exponentially modified Gaussian distribution 
% http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution

f = lambda * exp( ( sigma * lambda )^2 / 2 - ( x - mu ) * lambda ) .* normcdf( ( x - mu ) / sigma - sigma * lambda );
if sum( ~isfinite( f ) ) ~= 0 % the distribution is very close to a Gaussian, replace with a Gaussian
    f = 1 / ( sigma * sqrt( 2 * pi ) ) * exp( -( x - mu ).^2 / ( 2 * sigma^2 ) );
end