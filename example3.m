function [ dist, ld, dv ] = example3( pupil_diameter )
%
% demonstrates accommodation of the human eye by minimizing the retinal image
%
% Copyright: Yury Petrov, 2016
%

if nargin < 1
    pupil_diameter = 3; % mm, normal illumination
end

tic;
disp( 'Calculating...' );

nrays = 1000;                % execution time depends very little on the number of rays traced
nd = 30;                     % number of distances tested
dist = logspace( 2, 4, nd ); % tested pupil diameters, mm
%dist = linspace( 70, 1000, nd ); % tested pupil diameters, mm

ld = zeros( nd, 1 );
dv = zeros( nd, 1 );

% create an optical bench
bench = Bench;
eye = Eye( 10., pupil_diameter ); % eye with pupil diameter 10 mm
bench.append( eye );

for i = 1 : nd
    % create some rays
    if i < nd % diverging rays from a point source
        rays_in = Rays( nrays, 'source', [ -(dist( i ) - eye.Cornea1x) 0 0 ], [ 1 0 0 ], 1/0.82 * pupil_diameter / ( dist( i ) - eye.Cornea1x ), 'hexagonal','air' );
    else % the last point source is at infinity
        rays_in = Rays( nrays, 'collimated', [ -(dist( i ) - eye.Cornea1x) 0 0 ], [ 1 0 0 ], 1/0.82 * pupil_diameter, 'hexagonal','air' );
    end
    
    % find the tightest focus
    if i == 1
        lens_diameter = 10;
    else
        lens_diameter = ld( i - 1 ); % use the previous cycle value
    end
    options = optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'Diagnostics', 'off' );
    [ ld( i ), dv( i ) ] = fminunc( @focus, lens_diameter, options, bench, eye, rays_in );
end

figure, hold on;
plot( 0.1 * dist( 1 : end - 1 ), ld( 1 : end - 1 ), '-*' );
plot( 0.1 * dist( end ), ld( end ), '*' );
xlabel( 'Object distance, cm', 'FontSize', 18 );
ylabel( 'Lens diameter, mm', 'FontSize', 18 );
lbs = get( gca, 'XTickLabel' );
lbs{ end, : } = 'Inf ';
set( gca, 'XTickLabel', lbs );

figure, hold on;
plot( 0.1 * dist( 1 : end - 1 ), 1000 * dv( 1 : end - 1 ), '-o' );
plot( 0.1 * dist( end ), 1000 * dv( end ), 'o' );
xlabel( 'Object distance, cm', 'FontSize', 18 );
ylabel( 'Bundle focus, microns', 'FontSize', 18 );
ylim( [ 0 4 ] );
lbs = get( gca, 'XTickLabel' );
lbs{ end, : } = 'Inf ';
set( gca, 'XTickLabel', lbs );

[ mdv, mi ] = min( dv );
fprintf( 'Bundle tightests focus: %.5f mm\n', mdv );
fprintf( 'Optimal lens diameter: %.3f mm\n', ld( mi ) );
fprintf( 'Optimal distance: %.3f cm\n', 0.01 * dist( mi ) );

toc;


function dv = focus( ld, bench, eye, rays_in )
% trace the rays and find standard deviation of the image
eye.Lens( ld ); % change eye lens parameters
bench.elem{ 4 } = eye.elem{ 4 };
bench.elem{ 5 } = eye.elem{ 5 };
rays_through = bench.trace( rays_in );
[ ~, dv ] = rays_through( end ).stat; % record standard deviation of the bundle on the retina
%fprintf( '%.3f\t%.5f\n', ld, dv );
