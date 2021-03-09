function ellipse_draw( x0, cv, ax, angle, nvert, color )
% Draw an ellips defined by its center, half-axes, and angle
%
% INPUT:
% - x0 : 1 x 2 vector of the ellipse center
% - cv : 2 x 2 covariance matrix, or 1 x 3 vector of (1,1), (1,2), (2,2) elements of cv
% - ax : 1 x 2 vector of the ellipse principal half axes
% - ang : the angle of the ellipsoid rotation (longest half-axis off the
% x-axis, radians)
% 
% Copyright: Yury Petrov, 2016
%

if nargin < 6
    color = [ 0 0 1 ];
end
if nargin < 5
    nvert = 40;
end
if nargin < 4
    angle = [];
end
if nargin < 3
    ax = [];
end

nvert = ceil( nvert / 2 );

if size( cv, 1 ) == 2 && size( cv, 2 ) == 2
    cv = [ cv( 1, 1 ) cv( 1, 2 ) cv( 2, 2 ) ];
end

if isempty( ax ) || isempty( angle )
    d = sqrt( ( cv(1) - cv(3) )^2 / 4 + cv(2)^2 ); % discriminant
    t = ( cv(1) + cv(3) ) / 2;
    ax = 2 * sqrt( [ t + d; t - d ] ); % eigenvalues: [ max min ]
    angle = atan2( d - ( cv(1) - cv(3) ) / 2, cv(2) ); % positions the largest eigenvector along the x-axis
end

x = linspace( -ax(1), ax(1), nvert )';
y = ax(2) * sqrt( 1 - ( x / ax(1) ).^2 );
x = [ x;  flipud( x ) ];
y = [ y; -flipud( y ) ];
r = [ x y ];

r = r * [ cos( angle ), sin( angle ); -sin( angle ), cos( angle ) ];
r = r + repmat( x0, size( r, 1 ), 1 );
patch( r( :, 1 ), r( :, 2 ), color ); 