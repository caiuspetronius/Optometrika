function example14()
%
% Test astigmatic lens surfaces. The same as example1, but with an astigmatic front
% surface of the lens. The z (vertical dimension) curvature is changed here
% to produce vertical defocus.
%
% Copyright: Yury Petrov, 2017
%

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

% front lens surface
% parabolic surface with two different radii of curvature along the y and z
% dimensions (astigmatism). In this case the conic constant k applies to
% the y dimension, k value along the z dimension will be determined from the
% astigmatism Ry / Rz
lens1 = Lens( [ 40 0 0 ], 58, [ 40 45 ], -1, { 'air' 'bk7' } ); 
% back lens surface
lens2 = Lens( [ 60 0 0 ], 58, -70, -3, { 'bk7' 'air' } ); % concave hyperbolic surface
bench.append( { lens1, lens2 } );

% screen
screen = Screen( [ 102.76 0 0 ], 10, 10, 512, 512 ); % screen position based on the tightest convergence focus in example1
bench.append( screen );

%bench.rotate( [ 0 0 1 ], 0.15 );

% create some rays
nrays = 100;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 50, 'hexagonal','air' );

rays_through = bench.trace( rays_in );

% draw rays for the original focus (without the astigmatism)
rays_through = bench.trace( rays_in );    % repeat to get the min spread rays
bench.draw( rays_through, 'lines' );  % display everything, the other draw option is 'lines'

% get the screen image in high resolution
nrays = 10000;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 50, 'hexagonal','air' );
bench.trace( rays_in );
figure( 'Name', 'Image on the screen', 'NumberTitle', 'Off' );
imshow( screen.image, [] );

end

