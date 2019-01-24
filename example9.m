function example9()
%
% test cone mirrors
%
% Copyright: Yury Petrov, 2016
%

fov = 160; % FOV, degrees
fov = fov / 180 * pi;
D = [ 62 50 33 ]; % mirror diameters
LD = 47; % lens diameter

d1 = ( D(1) - D(2) ) / 2;
d2 = ( D(2) - D(3) ) / 2;
d3 = ( D(1) - D(3) ) / 2;
c = tan( fov / 4 );
r = roots( [ d3/2, c * ( d3/2 - d2 ), d3/2, c * ( d1 - d3/2 ) ] );
cang1 = -atan( r( imag( r ) == 0 ) ); %-fov/6;
cang2 = cang1 - fov/4;
h1 = -( D(1) - D(2) ) / 2 / tan( cang1 ); % depth of the first mirror
h2 = -( D(2) - D(3) ) / 2 / tan( cang2 ); % depth of the second mirror

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

block = Aperture( [ 30 0 0 ], [ 0 LD ] );
bench.append( block );

% first cone mirror surface
mirror1 = ConeLens( [ 30 0 0 ], D(1), h1, cang1, { 'air' 'mirror' } );
bench.append( mirror1 );

% second cone mirror surface
mirror2 = ConeLens( [ 30 + h1 0 0 ], D(2), h2, cang2, { 'air' 'mirror' } );
bench.append( mirror2 );

% screen
screen = Screen( [ 30 + h1 + h2 0 0 ], D(1), D(1), 256, 256 );
bench.append( screen );

% create collimated rays
nrays = 200;
rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], D(1) - 0.01, 'linear' );
y = rays_in.r( :, 2 );
z = rays_in.r( :, 3 );
rays_in.color( abs( z ) < 5 & y > LD/2, 2 ) = 0; % 
rays_in.color( abs( z ) < 5 & y > LD/2 & y <= LD/2 + ( D(1) - LD )/4, 1 ) = 1; % mark some rays red
rays_in.color( abs( z ) < 5 & y > LD/2 + ( D(1) - LD )/4, 3 ) = 1; % mark some rays blue

rays_through = bench.trace( rays_in, 0 ); % trace the rays, enable tracing rays that miss some bench elements by setting the second input parameter to 0

% draw bench elements and draw rays as arrows
bench.draw( rays_through, 'lines' );  % display everything, the other draw option is 'lines'

end
