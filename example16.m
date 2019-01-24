function example16()
%
% Export STL files of various lenses
%
% Copyright: Yury Petrov, 2018
%

nAngles = 100; % angular resolution of the STL export
flangeHeight = 2;

% polynomially aspheric lens
lens1 = AsphericLens( [ 0 0 0 ], 31, 80, 0, [ 0 -1e-05 3e-07 ], { 'air' 'pmma' } ); % front lens surface
lens2 = AsphericLens( [ 0 0 0 ], 31, -10, -2, [ 0 -2e-05 2e-07 ], { 'pmma' 'air' } ); % back lens surface
export_stl( lens1, lens2, 13, nAngles, [ flangeHeight flangeHeight ], 'aspheric.stl' );

% astigmatic ordinary lens
lens1 = Lens( [ 0 0 0 ], 58, [ 40 45 ], -1, { 'air' 'bk7' } ); % front lens surface
lens2 = Lens( [ 0 0 0 ], 58, -70, -3, { 'bk7' 'air' } ); % back lens surface
export_stl( lens1, lens2, 20, nAngles, [ flangeHeight flangeHeight ], 'astigmatic.stl' );

% general (cosine) lens
lens1 = GeneralLens( [ 0 0 0 ], 58, 'coslens', { 'air' 'bk7' }, 10, 116 ); % front lens surface
lens2 = Lens( [ 0 0 0 ], 58, -50, -3, { 'bk7' 'air' } ); % back lens surface
export_stl( lens1, lens2, 20, nAngles, [ flangeHeight flangeHeight ], 'Cosine.stl' );

% Fresnel lens
% front lens surface
R = 35; % lens radius of curvature
k = -1; % lens aspheric constant
D = 58; % lens diameter
nrings = 15; % number of the Fresnel rings
rads_outer = linspace( 5, D / 2, nrings ); % outer radii, the central ring's outer radius is 5 mm
rads_inner = [ 0 rads_outer( 1 : end - 1 ) ]; % the corresponding inner radii
d = sqrt( 1 - ( 1 + k ) * ( rads_inner / R ).^2 );
angs = pi/2 - atan( rads_inner ./ ( R * d ) ); % cone half-angle at the inner ring radius with respect to the x-axis direction
angs(1) = k; % replace the angle at the center (normally, pi/2) with the central part's conic constant k
pars = [ angs; ...
         repmat( pi - pi/8, 1, length( angs ) ); ... % wall half-angle at the outer ring radius, here less vertical by pi/8 (optional, defaults to pi)
         R ./ rads_outer ]; % ring's radial radius of curvature in units of the ring's outer radius (optional, defaults to Inf, i.e. to a flat conical surface)
sags = zeros( 1, length( rads_inner ) ); % collapse the lens along its axis to a Fresnel structure
lens1 = FresnelLens( [ 0 0 0 ], rads_outer, sags, pars, { 'air' 'pmma' } ); % Fresnel surface piece-wise identical to the surface of the original lens
% back lens surface
R2 = -115;
lens2 = Lens( [ 0 0 0 ], D, R2, 0, { 'pmma' 'air' } ); % concave hyperbolic surface
export_stl( lens1, lens2, 8, nAngles, [ flangeHeight flangeHeight ], 'Fresnel.stl' );
