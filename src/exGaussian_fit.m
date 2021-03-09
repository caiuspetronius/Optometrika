function [ mu, sigma, lambda ] = exGaussian_fit( x )
% EXGAUSSIAN_FIT fits exponentially modified Gaussian distribution to x
% http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
% 
% Copyright: Yury Petrov, 2016
%

m = mean( x );
s = std( x );
g = skewness( x );
sn = 1;
if g < 0
    g = -g;
    sn = -1;
end

t = ( g / 2 )^(1/3);
mu = m - s * t;
sigma = s * sqrt( 1 - t^2 );
lambda = sn / ( s * t ); % save sign into the sign of lambda