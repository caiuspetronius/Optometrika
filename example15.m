function example15()
%
% Simulate a hexagonal array of spherical microlenses
%
% Copyright: Yury Petrov, 2017
%

%Add src to PATH
startup;

R = 1;  % spherical lens radius
nLenses = 19;
D = 10 * R; % substrate plate diameter
h = 0.5;      % substrate plate height

% create a container for optical elements (Bench class)
bench = Bench;

% get hexagonal locations from the Rays function with the 'hexagonal' pattern
locs = Rays( nLenses, 'collimated', [ 10 0 0 ], [ 1 0 0 ], 7 * R, 'hexagonal','air' );

for i = 1 : length( locs.r )
    lens = Lens( locs.r( i, : ), 2 * R - 1e-6, R, 0, { 'air' 'pmma' } ); % a half-sphere
    bench.append( lens );
end

% add the screen-side substrate plate for the microlenses
plane = Plane( [ 10 + R + h 0 0 ], D + R, { 'pmma' 'air' } );
bench.append( plane );

% screen
screen = Screen( [ 12.3 0 0 ], D, D, 500, 500 );
bench.append( screen );

% create some rays
nrays = 300;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 10 * R, 'hexagonal','air' );
rays_out = bench.trace( rays_in, 0 );  % must use 0 as the last argument to consider rays missing a given microlens as candidates for hitting another one
bench.draw( rays_out, 'rays' );

% get a high-res image on the screen
rays_in = Rays( 100 * nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 10 * R, 'hexagonal','air' );
rays_out = bench.trace( rays_in, 0 );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( screen.image, [] );

