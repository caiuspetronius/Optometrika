function example8()
%
% test a lens with polynomial aspheric terms
%
% Copyright: Yury Petrov, 2016
%

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% front lens surface
lens1 = AsphericLens( [ 40 0 0 ], 31, 80, 0, [ 0 -1e-05 3e-07 ], { 'air' 'pmma' } );
bench.append( lens1 );

% back lens surface
lens2 = AsphericLens( [ 52.5 0 0 ], 31, -10, -2, [ 0 -2e-05 2e-07 ], { 'pmma' 'air' } ); % aspheric surface
bench.append( lens2 );

% screen
screen = Screen( [ 70 0 0 ], 3, 3, 256, 256 );
bench.append( screen );

% create collimated rays with some slant
nrays = 500;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 30, 'hexagonal','air' );

tic;
fprintf( 'Tracing rays... ' );
rays_through = bench.trace( rays_in );    % repeat to get the min spread rays

% draw bench elements and draw rays as arrows
bench.draw( rays_through );  % display everything, the other draw option is 'lines'

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 30, 'hexagonal','air' );
bench.trace( rays_in );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( screen.image, [] );

toc;

draw_lens_engineering( [ 40 52.5 ], [ 31 31 ], [ 80 -10 ], [ 0 -2 ], [ 0 -1e-05 3e-07; 0 -2e-05 2e-07 ]', 'pmma' ); 
fprintf( 'Engineering drawing saved in draw_lens_engineering.pdf\n' );
end
