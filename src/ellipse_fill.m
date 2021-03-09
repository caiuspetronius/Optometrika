function f = ellipse_fill( X, Y, x0, cv, ax, angle, hard )
% Fill in an ellips defined by its center, half-axes, and angle
%
% INPUT:
% - x0 : 1 x 2 vector of the ellipse center
% - cv : 2 x 2 covariance matrix or 1 x 3 vector of the upper triangular matrix of cv
% - ax : 1 x 2 vector of the ellipse principal half axes
% - ang : the angle of the ellipsoid rotation (longest half-axis off the x-axis, radians)
% - hard : a float indicating the ellipse edge hardness, 1 - makes ellispe similar to a 2D Gaussian, > 20 - a hard-edged ellipse
%
% OUTPUT: 
% - f : n x 1 matrix of values, 1 inside the ellipse, 0 - outside
% 
% Copyright: Yury Petrov, 2016
%

if nargin < 7
    hard = 100; % hard edged ellipse by default
end
if nargin < 6
    angle = [];
end
if nargin < 5
    ax = [];
end

if size( cv, 1 ) == 2 && size( cv, 2 ) == 2
    cv = [ cv( 1, 1 ) cv( 1, 2 ) cv( 2, 2 ) ];
end

if isempty( ax ) || isempty( angle )
    d = sqrt( ( cv(1) - cv(3) )^2 / 4 + cv(2)^2 ); % discriminant
    t = ( cv(1) + cv(3) ) / 2;
    ax = 2 * sqrt( [ t + d; t - d ] ); % eigenvalues: [ max min ]
    ax( ax < 0.5 ) = 0.5; % fill in ellipse dimensions less than a pixel. This is important to prevent 'tearing' of lens images
    angle = atan2( d - ( cv(1) - cv(3) ) / 2, cv(2) ); % positions the largest eigenvector along the x-axis
end

% rotation matrix
M = [ cos( -angle ), sin( -angle ); -sin( -angle ), cos( -angle ) ];

R = [ X(:) Y(:) ];
Rt = R - repmat( x0, size( R, 1 ), 1 ); % move to the origin
Rt = Rt * M; % rotate to align long half-axis with x-axis

f = zeros( size( Rt, 1 ), 1 );
r2 = ( Rt( :, 1 ) / ax(1) ).^2 + ( Rt( :, 2 ) / ax(2) ).^2;
if hard > 20
    ind = r2 < 2;
    f( ind ) = 1 ./ r2( ind ).^4; % antialiasing
    f( r2 <= 1 ) = 1; % hard-edged ellipse
else  % soft-edged ellipse
    ind = r2 < 1.5;
    f( ind ) = 0.5 * erfc( hard * ( sqrt( r2( ind ) ) - 1 ) );
end

if sum( f ) == 0 % the ellipse is smaller than a pixel
    [ ~, mi ] = min( r2 );
    f( mi ) = 1; % make the pixel closest to the center of the ellipse equal 1
end

f = reshape( f, size( X ) );