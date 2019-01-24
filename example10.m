function example10()
%
% test cylinder and cone surfaces with double refraction
%
% Copyright: Yury Petrov, 2016
%

D = 40; % cylinder diameter
h = 20; % cylinder height

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% cylindrical mirror surface with an elliptic cylinder
% the elliptic cylinder is defined by [ Dy Dz ] pair in the function below
mirror = CylinderLens( [ 20 0 0 ], [ D, 1.5*D ], h, { 'air', 'mirror' } );
bench.append( mirror );

% NOTE THAT THE CYLIDRICAL LENS IS DOUBLED HERE, BECAUSE RAYS
% ENCOUNTER IT TWICE!!!

% cylindrical lens front surface
lens = CylinderLens( [ 70 0 2*h ], D, 4*h, { 'air', 'pc' } );
lens.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens );
% cylindrical lens back surface
lens = CylinderLens( [ 70 0 2*h ], D, 4*h, { 'air' 'pc' } );
lens.rotate( [ 0 1 0 ], pi/2 );
bench.append( lens );

% conical lens front surface
% the elliptic cone is defined by [ Dy Dz ] pair in the function below
lens = ConeLens( [ 110 -55 0 ], [ 0.6*D D ], 108, 0.08, { 'air', 'pmma' } );
%lens = CylinderLens( [ 110 -55 0 ], [ 0.6*D D ], 108, { 'air', 'pmma' } );
lens.rotate( [ 0 0 1 ], pi/2 );
bench.append( lens );
% conical lens back surface
lens = ConeLens( [ 110 -55 0 ], [ 0.6*D D ], 108, 0.08, { 'air' 'pmma' } );
%lens = CylinderLens( [ 110 -55 0 ], [ 0.6*D D ], 108, { 'air', 'pmma' } );
lens.rotate( [ 0 0 1 ], pi/2 );
bench.append( lens );

% screen
screen = Screen( [ 150 0 0 ], 3 * D, 3 * D, 256, 256 );
bench.append( screen );

% create divergent rays
nrays = 100;
rays_in = Rays( nrays, 'source', [ 0 0 0 ], [ 1 0 0 ], 2., 'hexagonal' );

rays_through = bench.trace( rays_in, 0 ); % trace the rays, enable tracing rays that miss some bench elements by setting the second input parameter to 0

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'rays' );  % display everything, the other draw option is 'lines'

end
