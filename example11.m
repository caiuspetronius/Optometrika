function example11()
%
% Demonstrates ray tracing for rays originating inside the human eye
%
% Copyright: Yury Petrov, 2016
%

%Add src to PATH
startup;

bench = Bench;

eye = Eye;
eye.rotate( [ 0 0 1 ], pi ); % make the eye face along the positive x-axis direction
eye.remove( 1 ); % remove the retina surface, otherwise it will block the light even though it is physically behind the source
bench.append( eye );

screen = Screen( [ 20 0 0 ], 20, 20, 500, 500 );
bench.append( screen );

rays_in = Rays( 100, 'source', [ 0 0 0 ], [ 1 0 0 ], 0.5, 'hexagonal', 'vitreous' ); % make sure that the rays start with the right refractive indices, 'vitreous' in this case

rays_out = bench.trace( rays_in );
bench.draw( rays_out, 'rays' ); % draw results using dotted lines, the asterisks mark where rays intersect surfaces.