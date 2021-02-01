function example4()
%
% test a ring lens with the cosine surface profile defined in coslens.m
%
% Copyright: Yury Petrov, 2016
%


%Add src to PATH
startup;

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% aperture
% aper = Aperture( [ 40 0 0 ], [ 55 80 ] );
aper = Aperture( [ 40 0 0 ], [ 30 50 80 80 ] ); % rectangular aperture
bench.append( aper );

% front lens surface
lens1 = GeneralLens( [ 40 0 0 ], [ 20 58 ], 'coslens', { 'air' 'bk7' }, 10, 116 ); % cosine ring surface
% back lens surface
lens2 = Lens( [ 58 0 0 ], [ 20 58 ], -50, -3, { 'bk7' 'air' } ); % concave hyperbolic ring surface
bench.append( { lens1, lens2 } );

% prism front surface
prism1 = Plane( [ 65 0 0 ], [ 10, 30 ], { 'air' 'acrylic' } );
prism1.rotate( [ 0 0 1 ], 0.3 ); % rotate the front prism plane by .3 radians about z-axis
% prism back surface
prism2 = Plane( [ 75 0 0 ], [ 10, 30 ], { 'acrylic' 'air' } );
prism2.rotate( [ 0 1 0 ], 0.2 ); % rotate the back prism plane by .2 radians about y-axis
bench.append( { prism1, prism2 } );

% screen
screen = Screen( [ 93 0 0 ], 20, 15, 128 * 20/15, 128 );
bench.append( screen );

% create collimated rays with some slant
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 -0.1 0 ], 68, 'hexagonal','air' );

tic;
fprintf( 'Tracing rays... ' );
rays_through = bench.trace( rays_in );

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'arrows' );  % display everything, the other draw option is 'lines'

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 -0.1 0 ], 58, 'hexagonal','air' );
bench.trace( rays_in, 0 );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( kron( screen.image, ones( 3 ) ), [] );

toc;
end
