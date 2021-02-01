function [ mu, sigma, lambda, ang ] = exGaussian_fit2D( x )
% EXGAUSSIAN_FIT2D fits exponentially modified Gaussian distribution to x
% in 2D
% http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
% 
% Copyright: Yury Petrov, 2016
%

cv = cov( x );
d = sqrt( ( cv(1,1) - cv(2,2) )^2 / 4 + cv(1,2)^2 ); % discriminant
t = ( cv(1,1) - cv(2,2) ) / 2;
V = [ cv(1,2), -cv(1,2); d - t, d + t ]; % eigenvectors
V( :, 1 ) = V( :, 1 ) / norm( V( :, 1 ) ); % normalize
V( :, 2 ) = V( :, 2 ) / norm( V( :, 2 ) );
xt = x * V;  % rotate data into eigevectors RF, the largest eigenvector is along the x-axis
ang = atan2( V( 2, 1 ), V( 1, 1 ) ); % positions the largest eigenvector along the x-axis

% fit each dimension separately
[ mu( 1 ), sigma( 1 ), lambda( 1 ) ] = exGaussian_fit( xt( :, 1 ) );
[ mu( 2 ), sigma( 2 ), lambda( 2 ) ] = exGaussian_fit( xt( :, 2 ) );

if lambda( 1 ) < 0
    ang = ang - pi;
    xt = -xt; % rotate by Pi
    [ mu( 1 ), sigma( 1 ), lambda( 1 ) ] = exGaussian_fit( xt( :, 1 ) );
    [ mu( 2 ), sigma( 2 ), lambda( 2 ) ] = exGaussian_fit( xt( :, 2 ) );
    
    if lambda( 2 ) < 0 % the sign of lambda remains and indicates the flip
        xt( :, 2 ) = -xt( :, 2 ); % reverse y sign
        mu( 2 ) = exGaussian_fit( xt( :, 2 ) );
    end
end
