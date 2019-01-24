function example12()
%
% Draw a lens and determine its front surface, back surface, and total
% height. Make an animated gif of the lens and an engineering drawing of
% the lens.
%
% Copyright: Yury Petrov, 2016
%

Df = 40; % front diameters
Db = 50; % back diameters
Rf = -50; % front radius of curvature
Rb = 200; % back radius of curvature
kf = -3; % front conic constant
kb = 0; % back conic constant
h = 14; % lens height at the edge

[ ht, hf, hb, V ] = lens_dims( Df, Db, Rf, Rb, kf, kb, h, [], [], 1 );
fprintf( 'Front surface height: %.3f\n', hf );
fprintf( 'Back  surface height: %.3f\n', hb );
fprintf( 'Total lens height in the center: %.3f\n', ht );
fprintf( 'Lens volume: %.3f\n', V );
fprintf( 'Animated gif saved in lens_dims.gif\n' );
f = draw_lens_engineering( [ 0, hf + hb + h ], [ Df Db ], [ Rf Rb ], [ kf kb ], [], 'zeonex', [], 'AR', 'flange', 2, 2, 'Make my lens really beautiful');
fprintf( 'Engineering drawing saved in draw_lens_engineering.pdf\n' );