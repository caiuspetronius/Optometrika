function [ x, yp ] = frenlens( y, Fr )
% FRENLENS defines a radially symmetric Fresnel profile.
% the first input argument is the profile radii, the
% second argument is a FresnelLens object.
%
% Copyright: Yury Petrov, 2017
%

x = [];
yp = [];
for i = 1 : Fr.ncones % loop over cones
    if i == 1
        if y(1) == 0
            radin = 0;
        else
            2 * Fr.rad( 1, 1 ) - Fr.rad( 2, 1 ); % take the inner radius assuming the same step
        end
    else
        radin = Fr.rad( i - 1, 1 );
    end
    
    if isempty( Fr.vx ) || Fr.R( i, 1 ) == 0 || isinf( Fr.k( i ) ) % cone surface
        z( 1, : ) = Fr.sag( i );
        z( 2, : ) = Fr.sag( i ) + ( Fr.rad( i, 2 ) - radin ) / tan( Fr.the( i, 1 ) );
    else % quadric surface
        r = y( y >= radin & y < Fr.rad( i, 2 ) );
        a = 1 + Fr.k( i );
        r2 = r.^2 / Fr.R( i )^2;
        if abs( a ) < 1e-10 % paraboloid, special case
            z = r2 * Fr.R( i, 1 ) / 2;
        else
            z = Fr.R( i ) / a * ( 1 - sqrt( 1 - a * r2 ) );
        end
        z = z + Fr.sag( i ) + Fr.vx( i ); % add sag and quadric vertex displacement
    end
    
    if abs( Fr.the( i, 1 ) - pi/2 ) > 1e-10 && ( sign( ( z( 2 ) - z( 1 ) ) / ( Fr.rad( i, 2 ) - radin ) ) ~= sign( tan( Fr.the( i, 1 ) ) ) )
        z = flipud( z );
    end
    x = [ x z ];
    yp = [ yp r ];
end
