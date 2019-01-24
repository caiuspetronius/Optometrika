function rinter = conic_intersection( r_in, e, surf )
%
% returns intersection vector with ellipsoid, paraboloid, or hyperboloid
% surface
%

x0 = r_in( :, 1 );
y0 = r_in( :, 2 );
z0 = r_in( :, 3 );
e1 = e( :, 1 );
e2 = e( :, 2 );
e3 = e( :, 3 );
k = surf.k;
a = 1 + k;
R = surf.R(1);

if a == 0 % paraboloid, special case
    A = e2.^2 + e3.^2;
    B = e1 * R - e2 .* y0 - e3 .* z0;
    D = B.^2 - A .* ( -2 * R * x0 + y0.^2 + z0.^2 );
    D( D < 0 ) = Inf; %  mark no intersection as Inf, I is nulled below anyway
    
    d01 = ( y0.^2 + z0.^2 ) / ( 2 * R ) - x0; % distance to the intersection for the ray || to the paraboloid
    d02 = d01;
else
    A = e1 .* ( a * e1.^2 + e2.^2 + e3.^2 );
    B = e1.^2 * R - a * e1.^2 .* x0 - e1 .* e2 .* y0 - e1 .* e3 .* z0;
    D = e1.^2 .* ( e3.^2 .* ( 2 * R * x0 - a * x0.^2 - y0.^2 ) + 2 * e2 .* e3 .* y0 .* z0 - ...
        2 * e1 .* ( R - a * x0 ) .* ( e2 .* y0 + e3 .* z0 ) + e2.^2 .* ( 2 * R * x0 - a * x0.^2 - z0.^2 ) + ...
        e1.^2 .* ( R^2 - a * ( y0.^2 + z0.^2 ) ) );
    D( D < 0 ) = Inf; %  mark no intersection as Inf, I is nulled below anyway
    
    A0 = 2 * ( 2 * a * e2 .* e3 .* y0 .* z0 + e2.^2 .* ( R^2 - a * z0.^2 ) + e3.^2 .* ( R^2 - a * y0.^2 ) );
    B0 = a * ( 1 - a ) * ( a * ( e1 .* x0 - e2 .* y0 - e3 .* z0 ) - R * e1 ) .* ( y0.^2 + z0.^2 );
    d01 =  B0 ./ A0;   % distance to the intersection for the ray || to the hyperboloid sides
    d02 = -B0 ./ A0;
end
d1 = ( B + sqrt( D ) ) ./ A;
d2 = ( B - sqrt( D ) ) ./ A;

% it is necessary to eliminate infinities before the following logical operation
d1(  ~isfinite( d1 )  ) = 0;
d2(  ~isfinite( d2 )  ) = 0;
d01( ~isfinite( d01 ) ) = 0;
d02( ~isfinite( d02 ) ) = 0;

d( :, 1 ) = ( abs( A ) <= eps ) .* d01 + ... % ray parallel to the paraboloid / hyperboloid, special case
    ( abs( A ) >  eps ) .* d1;
d( :, 2 ) = ( abs( A ) <= eps ) .* d02 + ... % ray parallel to the paraboloid / hyperboloid, special case
    ( abs( A ) >  eps ) .* d2;

% find the shortest positive distance to the (two) intersections along the ray
if a * R < 0 % hyperboloid with positive R, or ellipsoid with negative R: we want to consider only the top branch
    ind = x0 - R / a * ( 1 + sqrt( 1 - a * ( y0.^2 + z0.^2 ) / R^2 ) ) < 0; % sources below the bottom branch
    d( ind, 1 ) = NaN; %realmax; % disregard intersections with the lower branch
elseif a * R > 0 % hyperboloid with negative R, or ellipsoid with positive R: we want to consider only the bottom branch
    ind = x0 - R / a * ( 1 - sqrt( 1 - a * ( y0.^2 + z0.^2 ) / R^2 ) ) > 0; % sources above the bottom branch
    d( ind, 1 ) = NaN; %realmax; % disregard intersections with the lower branch
end
d( d <= 1e-12 ) = NaN; %realmax; % intensities for these rays (non-intersecting the surface) will be set to 0 anyway

[ ~, ii ] = min( d, [], 2 );
ind = sub2ind( size( d ), ( 1:size( d, 1 ) )', ii ); % linear index of the min value
d = abs( d( ind ) );

% form the intersection vector
rinter = r_in + repmat( d, 1, 3 ) .* e;
                        
