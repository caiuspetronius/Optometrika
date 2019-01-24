function example2()
%
% demonstrates the Optometrika's optical model of the human eye
%
% Copyright: Yury Petrov, 2016
%

% demonstrate resolution dependence on the pupil diameter
tic;
disp( 'Calculating...' );

% add some optical elements
nrays = 1000;
nd = 20;                % number of pupil diameters tested
pd = linspace( 2, 8, nd ); % tested pupil diameters, mm

dv = zeros( nd, 1 );
for i = 1 : nd
    % create an optical bench
    bench = Bench;
    pupil_diameter = pd( i );
    % eye = Eye( [], pupil_diameter ); % eye accommodating at infinity
    eye = Eye( 9.4, pupil_diameter ); % eye accommodating at infinity
    %eye.rotate( [ 0 0 1 ], 1=30 * pi / 180 );
    bench.append( eye );
    
    % create some rays
    rays_in = Rays( nrays, 'collimated', [ -20 0 0 ], [ 1 0 0 ], 1.2 * pupil_diameter, 'hexagonal' );
    
    % trace the rays
    rays_through = bench.trace( rays_in );
    [ ~, dv( i ) ] = rays_through( end ).stat; % record standard deviation of the bundle on the retina
end

[ mdv, mi ] = min( dv );
fprintf( 'Bundle tightests focus: %.5f mm\n', mdv );
fprintf( 'Optimal pupil diameter: %.5f mm\n', pd( mi ) );

toc;
tic;
disp( 'Drawing...' );

figure( 'Name', 'Eye focusing vs. Pupil diameter', 'NumberTitle', 'Off' );
plot( pd, dv, '*-' );
xlabel( 'Pupil diameter, mm', 'FontSize', 16 );
ylabel( 'Bundle focus (standard deviation), mm', 'FontSize', 16 );

% draw the result for the optimal pupil diameter
bench = Bench;
eye = Eye( [], pd( mi ) );
% eye.rotate( [ 0 0 1 ], 20 * pi / 180 );
bench.append( eye );
nrays = 40; 
rays_in = Rays( nrays, 'collimated', [ -20 0 0 ], [ 1 0 0 ], 1.2 * pd( mi ), 'hexagonal' );
rays_in.color( ceil( rays_in.cnt / 2 ), : ) = [ 1 0 0 ]; % draw the central ray in red
rays_through = bench.trace( rays_in );
bench.draw( rays_through, 'lines' );

toc;

% accommodation (lens diameter) from infinity to 7D (14 cm)
tic;
disp( 'Calculating...' );

% add some optical elements
pupil_diameter = 3; % normal light level
nrays = 1000;
nd = 30;                % number of lens diameters tested
ld = linspace( 10.610, 9.5, nd ); % tested pupil diameters, mm
viewing_distance = 140; % mm

% create some rays
rays_in = Rays( nrays, 'source', [ -(viewing_distance + 13.3) 0 0 ], [ 1 0 0 ], 1.2 * pupil_diameter / (viewing_distance + 13.3), 'hexagonal' );
    
dv = zeros( nd, 1 );
for i = 1 : nd
    % create an optical bench
    bench = Bench;
    lens_diameter = ld( i );
    eye = Eye( lens_diameter, pupil_diameter ); % eye accommodating at infinity
    %eye.eye_lens_vol( lens_diameter )
    bench.append( eye );    
    % trace the rays
    rays_through = bench.trace( rays_in );
    [ ~, dv( i ) ] = rays_through( end ).stat; % record standard deviation of the bundle on the retina
end

[ mdv, mi ] = min( dv );
fprintf( 'Bundle tightests focus: %.5f mm\n', mdv );
fprintf( 'Optimal lens diameter: %.5f mm\n', ld( mi ) );

toc;
tic;
disp( 'Drawing...' );

figure( 'Name', 'Eye focusing vs. Lens diameter', 'NumberTitle', 'Off' );
plot( ld, dv, '-o' );
xlabel( 'Lens diameter, mm', 'FontSize', 16 );
ylabel( 'Bundle focus (standard deviation), mm', 'FontSize', 16 ); 
set( gca, 'XDir', 'reverse' );

% draw the result for the optimal pupil diameter
bench = Bench;
eye = Eye( ld( mi ), pupil_diameter );
bench.append( eye );
nrays = 40; 
rays_in = Rays( nrays, 'source', [ -(viewing_distance + 13.3) 0 0 ], [ 1 0 0 ], 1.2 * pupil_diameter / (viewing_distance + 13.3), 'hexagonal' );
rays_through = bench.trace( rays_in );
bench.draw( rays_through, 'lines' );

toc;

end

