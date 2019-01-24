function s = aspheric( pars, r )
%
% ASPHERIC describes a symmetric aspherical lens with sag given by
% s = r^2 / R / ( 1 + sqrt( 1 - ( k + 1 ) * r^2 / R^2 ) ) + a(1) * r^2 +
% a(2) * r^4 + ... + a(n) * r^2n
%
% INPUT:
%    pars: with the following parameters
%       pars(1:2) = R - two radii of curvature
%       pars(3) = k - conic constant
%       pars(4:end) = a - higher aspherics values
%    r - radii
%
% OUTPUT:
%    s - sag values corresponding to r values
%
% Copyright: Yury Petrov, 2016
%

R = pars( 1 ); % radius of curvature along y-axis
k = pars( 3 );
a = pars( 4:end );
r2 = r.^2;
s = ( r2 / R ) ./ ( 1 + sqrt( 1 - ( k + 1 ) * r2 ./ R^2 ) );
s( ~isreal( s ) ) = 1e+20; % prevent complex values
for i = 1 : length( a )
    s = s + a( i ) * r2.^i;
end