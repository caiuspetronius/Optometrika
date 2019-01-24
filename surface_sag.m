function s = surface_sag( D, R, k, a )
%
% calculate maximum sag of the surface defined by diameter D, radius R,
% and conic constant k, and, optionally, by aspheric polynomial terms a
%
% Copyright: Yury Petrov, 2016
%

if nargin < 4 || isempty( a ) || sum( a.^2 ) == 0 % no polynomial coefficients
    if k == -1 % paraboloid
        s = D^2 / ( 8 * R );
    else % spheroids, hyperboloids
        if isinf( R )
            s = 0;
        else
            x = R / ( 1 + k );
            det = 1 - ( 1 + k ) * D^2 / ( 4 * R^2 );
            if det < 0
                s = Inf;
            else
                s = sign( R ) * abs( x * ( 1 - sqrt( det ) ) );
            end
        end
    end
else
    if size( a, 1 ) > size( a, 2 )
        a = a';
    end
%     r = linspace( 0, D/2, 10000 ); % densly sample the lens radius
%     s = max( abs( aspheric( [ R k a ], r ) ) );
      s = aspheric( [ R k a ], D/2 );
end
