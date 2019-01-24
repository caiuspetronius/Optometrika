function [ x0, cv, ax, ang ] = ellipse_fit( x )
% fits an ellipse to cover a cloud of dots in 2D defined by x
%
% INPUT: 
% - x : n x 2 or n x 3 matrix of datapoints
%
% OUTPUT:
% - x0 : 1 x 2 vector of the ellipse center
% - ax : 1 x 2 vector of the ellipse principal half axes
% - ang : the angle of the ellipsoid rotation (longest half-axis off the
% x-axis, radians)
% 
% Copyright: Yury Petrov, 2016
%

ax = [];
ang = [];

if nargin < 2
    flag = 0;
end

x0 = mean( x );
cv = cov( x );  % get covariance matrix

if flag == 0 % use x moments
    if size( x, 2 ) == 2
        d = sqrt( ( cv(1,1) - cv(2,2) )^2 / 4 + cv(1,2)^2 ); % discriminant
        t = ( cv(1,1) + cv(2,2) ) / 2;
        ax = 2 * sqrt( [ t + d; t - d ] ); % eigenvalues: [ max min ]
        ang = atan2( d - ( cv(1,1) - cv(2,2) ) / 2, cv(1,2) ); % positions the largest eigenvector along the x-axis
    else % do PCA
        [ ~, S, V ] = svd( cv );
        ax = sqrt( diag( S ) );
        ax = ax( 1 : 2 ); % choose the two largest eigenvalues
        ang = atan2( V( 2, 1 ), V( 3, 1 ) );
    end
end
