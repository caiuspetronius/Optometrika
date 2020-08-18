function example7()
%
% test a Fresnel lens with quadric rings
%
% Copyright: Yury Petrov, 2017
%

% create a container for optical elements (Bench class)
bench = Bench;

condition = 'fresnel'; % change to 'smooth' to see the original lens or to 'fresnel-sim' to simulate the original lens with a globally smooth Fresnel surface



% front lens surface
% This is the simulated aspheric lens surface
R = 35; % lens radius of curvature
k = -1; % lens aspheric constant
D = 58; % lens diameter
nrings = 15; % number of the Fresnel rings

x1 = 30;
if strcmp( condition, 'smooth' )
    lens1 = Lens( [ x1 0 0 ], D, R, k, { 'air' 'pmma' } ); % smooth surface
else
    % This simulates the above surface with a Fresnel surface
    rads_outer = linspace( 5, D / 2, nrings ); % outer radii, the central ring's outer radius is 5 mm
    rads_inner = [ 0 rads_outer( 1 : end - 1 ) ]; % the corresponding inner radii
    d = sqrt( 1 - ( 1 + k ) * ( rads_inner / R ).^2 );
    angs = pi/2 - atan( rads_inner ./ ( R * d ) ); % cone half-angle at the inner ring radius with respect to the x-axis direction
    angs(1) = k; % replace the angle at the center (normally, pi/2) with the central part's conic constant k
     
%     % Fresnel surface made of cone segments
%     pars = angs; 
    
    % Fresnel surface made of quadric segments
    pars = [ angs; ... 
             repmat( pi - pi/8, 1, length( angs ) ); ... % wall half-angle at the outer ring radius, here less vertical by pi/8 (optional, defaults to pi)
             R ./ rads_outer ]; % ring's radial radius of curvature in units of the ring's outer radius (optional, defaults to Inf, i.e. to a flat conical surface)
    if strcmp( condition, 'fresnel' )
        sags = zeros( length( rads_inner ), 1 ); % collapse the lens along its axis to a Fresnel structure
        x1 = 42;
        lens1 = FresnelLens( [ x1 0 0 ], rads_outer, sags, pars, { 'air' 'pmma' } ); % Fresnel surface piece-wise identical to the surface of the original lens
    else
        if k == -1
            sags = rads_inner.^2 ./ ( 2 * R );
        else
            sags = R / ( 1 + k ) * ( 1 - d ); % keep the original lens's profile
        end
        lens1 = FresnelLens( [ x1 0 0 ], rads_outer, sags, pars, { 'air' 'pmma' } ); % Fresnel surface globally identical to the surface of the original lens
    end
end
bench.append( lens1 );

% back lens surface
R2 = -115;
x2 = 48;
lens2 = Lens( [ x2 0 0 ], D, R2, 0, { 'pmma' 'air' } ); % concave hyperbolic surface
bench.append( lens2 );

% screen
screen = Screen( [ 94.5 0 0 ], 5, 5, 512, 512 );
bench.append( screen );

tic;
fprintf( 'Tracing rays... ' );
nrays = 100;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 50, 'hexagonal' );
rays_through = bench.trace( rays_in );
toc;

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'rays' );  % display everything, the other draw option is 'lines'
figure, imshow( screen.image, [ 0 max( max( screen.image ) ) ] );

% draw the lens for manufacturing

if strcmp( condition, 'fresnel' )
    [ ~, pr ] = draw_lens_engineering( [ 0, x2 - x1 ], [ D D ], [ R R2 ], [ k 0 ], [], 'pmma', [], 'AR', 'Fresnel', lens1 );
else
    [ ~, pr ] = draw_lens_engineering( [ 0, x2 - x1 ], [ D D ], [ R R2 ], [ k 0 ], [], 'pmma', [], 'AR'); 
end

fprintf( 'Engineering drawing saved\n' );
save lens_profile.mat pr;


