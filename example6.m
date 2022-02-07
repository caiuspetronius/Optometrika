function example6()
%
% test planar and parabolic mirrors (a refractor telescope)
%
% Copyright: Yury Petrov, 2016
%

%Add src to PATH
startup;

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% back surface of the planar mirror
mirror1 = Plane( [ 45 0 0 ], 8, 4, { 'mirror' 'air' } ); % pay attention to the glass order here!
mirror1.rotate( [ 0 0 1 ],  -pi / 4 );
bench.append( mirror1 );

% parabolic mirror
mirror2 = Lens( [ 100 0 0 ], 40, -120, -1, { 'air' 'mirror' } ); % pay attention to the glass order here!
bench.append( mirror2 );

% front sufrace of the planar mirror
mirror3 = Plane( [ 45 0 0 ], 8, 4, { 'mirror' 'air' } ); % pay attention to the glass order here!
mirror3.rotate( [ 0 0 1 ],  -pi / 4 );
bench.append( mirror3 );

% screen
screen = Screen( [ 45 -5 0 ], 1, 1, 512, 512 );
screen.rotate( [ 0 0 1 ], -pi/2 );
bench.append( screen );

% create collimated rays
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 58, 'hexagonal','air' );

fprintf( 'Tracing rays...\n' );
rays_through = bench.trace( rays_in, 0 ); % the second parameter set to 0 enables ray tracing for rays missing some elmeents on the bench

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'arrows', [], [ 3 0 2 1 .1 ] );  % display everything, specify arrow length for each bench element

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 58, 'hexagonal','air' );
bench.trace( rays_in, 0 );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( screen.image, [] );

end
