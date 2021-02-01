function example5()
%
% test planar mirrors
%
% Copyright: Yury Petrov, 2016
%


%Add src to PATH
startup;


% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% planar mirror
mirror1 = Plane( [ 60 0 0 ], 40, 40, { 'air' 'mirror' } ); % pay attention to the glass order here, it defines the mirror orientation!
mirror1.rotate( [ 0 0 1 ], -pi / 4 );
bench.append( mirror1 );

% planar mirror
mirror2 = Plane( [ 60 50 0 ], 40, 40, { 'mirror' 'air' } ); % pay attention to the glass order here!
mirror2.rotate( [ 0 0 1 ], -pi / 4 );
bench.append( mirror2 );

% screen
screen = Screen( [ 100 50 0 ], 20, 15, 128 * 20/15, 128 );
bench.append( screen );

% create collimated rays
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 58, 'hexagonal','air' );

fprintf( 'Tracing rays...\n' );
rays_through = bench.trace( rays_in );

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'lines', [], 2 );  % display everything, scale arrow length 2x

end
