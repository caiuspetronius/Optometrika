function x = asphlens( y, z, args, flag )
% ASPHLENS defines a commonly used radially symmetric aspherical profile.
% the first two input arguments are coordinates in the lens plane, the
% third argument is a vector holding the aspheric parameters:
%       args(1:2) = Ry and Rz - radii of curvature
%       args(3) = k - conic constant along y
%       args(4:end) = a - higher aspherics values
% On flag == 0 the function 
% should return the lens height x for the given position ( y, z ). Otherwise,
% the function should return the lens normal at this position. By convention, 
% the normal should point along the x-axis, i.e. in the same general 
% direction as the traced ray.
%
% Copyright: Yury Petrov, 2017
%

R = args(1); % Ry
sc = args(1) / args(2); % astigmatism: Ry / Rz
k = args(3);
r = sqrt( y.^2 + z.^2 );
r2yz = ( y.^2 + z.^2 ) / R^2;
if flag == 0
    x = aspheric( args, sqrt( y.^2 + ( sc * z ).^2 ) ); % aspheric surface profile
else
    if k == -1 % parabola, special case
        c = 1 ./ sqrt( 1 + r2yz );
        s = sqrt( 1 - c.^2 );
    else
        s = sqrt( r2yz ) ./ sqrt( 1 - k * r2yz );
        c = sqrt( 1 - s.^2 );
    end
    tg = sign( R ) * s ./ c;
    for i = 4 : length( args ) % add higher aspherical derivatives to the conic derivative
        ii = i - 2;
        tg = tg + args( i ) * 2 * ii .* r.^( 2 * ii - 1 );
    end
    c = 1 ./ sqrt( 1 + tg.^2 ); % component along the ray direction, always positive
    s = -sign( tg ) .* sqrt( 1 - c.^2 ); % sign of the transverse component determined by the local convave/convex shape
    th = atan2( z, y ); % rotation angle to bring r into XZ plane
    x = [ c, s .* cos( th ), sc * s .* sin( th ) ]; % make normal sign positive wrt ray, scale z-component according to the one-form transformation
end

end
