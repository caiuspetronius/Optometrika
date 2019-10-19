function example13()
%
% test refraction through the lens edge and backward rays refraction using the sub-aperture Maksutov-Cassegrain
% telescope design
%
% Copyright: Yury Petrov, 2016
%

% create a container for optical elements (Bench class)
bench = Bench;
hD = 15; % diameter of the hole in the concave spherical mirror
D = 30;  % distance to the concave mirror
d = D - 28;  % distance to the convex mirror

% add optical elements in the order they are encountered by light rays

% back surface of the convex mirror
mirror1 = Lens( [ d 0 0 ], hD, -90, -1.1, { 'mirror' 'air' } ); % pay attention to the glass order here!
bench.append( mirror1 );

% spherical mirror
mirror2 = Lens( [ D 0 0 ], [ hD 40 ], -89.5, 0, { 'air' 'mirror' } ); % pay attention to the glass order here!
bench.append( mirror2 );

% meniscus lens on the way from concave to convex mirror
lens1 = Lens( [ d + 2 0 0 ], hD, 100, 0, { 'bk7' 'air' } ); % pay attention to the glass order here!
lens2 = Lens( [ d + 1 0 0 ], hD, 105, 0, { 'air' 'bk7' } ); % pay attention to the glass order here!
lens2sag = surface_sag( hD, 105, 0 );
cylin = CylinderLens( [ d + 1 + lens2sag 0 0 ], hD, 1, { 'bk7' 'air' } ); % cylindrical lens equator surface
bench.append( lens1 );
bench.append( cylin );
bench.append( lens2 );

% front surface of the convex mirror
mirror3 = Lens( [ d 0 0 ], hD, -90, -1.1, { 'mirror' 'air' } ); % pay attention to the glass order here!
bench.append( mirror3 );

% meniscus lens on the way from convex mirror to the screen
%
% FOR BOTH FORWARD (ALONG POSITIVE X-AXIS) AND BACKWARD (ALONG NEGATIVE X-AXIS) RAYS ORDER INTERFACE MATERIAL AS IF FOR FORWARD RAY
% DIRECTION!!! HENCE, CAN REUSE THE LENS SURFACES.
%
bench.append( lens2 ); % pay attention to the order of surfaces here
bench.append( lens1 );

% screen
screen = Screen( [ 29.34 0 0 ], 0.5, 0.5, 512, 512 );
bench.append( screen );

% create collimated rays
nrays = 100;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 40, 'hexagonal' );

fprintf( 'Tracing rays...\n' );
rays_through = bench.trace( rays_in, 0 ); % the second parameter set to 0 enables ray tracing for rays missing some elmeents on the bench

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'lines', [], [ 3 0 2 1 .1 ] );  % display everything, specify arrow length for each bench element
fp = rays_through( end ).focal_point;
fprintf( 'The focal point of the system at: %.3f\n', fp(1) );

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 40, 'hexagonal' );
bench.trace( rays_in, 0 );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( screen.image, [] );

end
