function [ r1, r2, t1, t2 ] = quadric_intersection( a2, b2, c2, s, r0, n )
%
% Find intersection of a line defined by the origin vector r0 and direction
% vector n with a quadric surface defined by a^2, b^2, c^2 (signs
% included) and the sign of the constant 1 (s), i.e.  with the equation
% x^2/a^2 + y^2/b^2 + z^2/c^2 = s
% https://en.wikipedia.org/wiki/Quadric
% a2 = 0 and s = 0 defines a paraboloid (elliptical or hyperbolic)
% a2 = 0 and s = 1 defines a cylinder (elliptical, hyperbolic, or parabolic)
%
% INPUT: 
% - a2, b2, c2 : quadric semiaxes. Sign of b2 should be positive
% - s : the constant of the quadric as given by x^2/a^2 + y^2/b^2 + z^2/c^2 = s
% - r0 : rays origins, n x 3 matrix
% - n : ray directions, n x 3 matrix
%
% OUTPUT:
% - r1, r2 : n x 3 intersection matrices
% - t1, t2 : n distances to the intersections
%
% Copyright: Yury Petrov, 2017
%

if b2 <= 0
    error( 'b2 must be > 0!' );
end

x0 = r0( :, 1 );
y0 = r0( :, 2 );
z0 = r0( :, 3 );
x02 = x0.^2;
y02 = y0.^2;
z02 = z0.^2;

j = n( :, 1 );
k = n( :, 2 );
l = n( :, 3 );
j2 = j.^2;
k2 = k.^2;
l2 = l.^2;

if a2 == 0
    if s == 1 % cylinder
        a = k2 ./ b2 + l2 ./ c2;
        b = 2 * ( y0 .* k ./ b2 + z0 .* l ./ c2 );
        c = y02 ./ b2 + z02 ./ c2 - 1;
    elseif s == 0 % paraboloid
        a = k2 ./ b2 + l2 ./ c2;
        b = 2 * ( y0 .* k ./ b2 + z0 .* l ./ c2 - j/2 );
        c = y02 ./ b2 + z02 ./ c2 - x0;
    end
else % ellipsoids, hyperboloids, cones
    a = j2 ./ a2 + k2 ./ b2 + l2 ./ c2;
    b = 2 * ( x0 .* j ./ a2 + y0 .* k ./ b2 + z0 .* l ./ c2 );
    c = x02 ./ a2 + y02 ./ b2 + z02 ./ c2 - s;
end

D = b.^2 - 4 * a .* c;
D( D < -1e-12 ) = Inf; % suppress negative values to avoid imaginary values
D( D < 0 ) = 0; % remove round off negative Ds
t1 = 0.5 * ( -b - sqrt( D ) ) ./ a;
t2 = 0.5 * ( -b + sqrt( D ) ) ./ a;

dims = size(n);
r1 = r0 + n .* repmat(t1,1,dims(2));
r2 = r0 + n .* repmat(t2,1,dims(2));

