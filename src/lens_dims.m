function [ ht, hf, hb, V ] = lens_dims( Df, Db, Rf, Rb, kf, kb, h, af, ab, disp_fl )
% LENS_DIMS calculates a lens thickness and volume given its surface
% parameters, draws the lens and save the animation in lens_dims.gif
% Volume is calculated correctly only for lenses without astigmatism.
%   
% INPUT:
%    Df, Db - lens diameter(s), front and back
%    Rf - radius of the tangent sphere for the front surface
%    Rb - radius of the tangent sphere for the back surface
%    kf - asphreicity coefficient for the front surface
%    kb - asphreicity coefficient for the back surface
%    h  - lens height at the equator (conical part)
%    af - aspheric (polynomial) terms for the front surface
%    ab - aspheric (polynomial) terms for the bakc surface
%    disp_fl - (1, default) display the lens and save its animation in lens_dims.gif, (0) do not
%    dispaly the lens
%
% OUTPUT:
%    ht - total lens thickness
%    hf - front part thickness
%    hb - back part thickness
%    V  - lens volume
%
% Copyright: Yury Petrov, 2016
%

if nargin < 10
    disp_fl = 1;
end
if nargin < 9
    ab = [];
end
if nargin < 8
    af = [];
end
if nargin < 7
    h = 0;
end

if length( Df ) == 1
    Df = [ 0 Df ];
end
if length( Db ) == 1
    Db = [ 0 Db ];
end


% get volumes and heights at the center
nr = 10000;
r = linspace( Df(1)/2, Df(2)/2, nr );
dr = ( Df(2) - Df(1) ) / 2 / nr;
sf = aspheric( [ Rf(1) Rf(1) kf af ], r );
hf = sf( end );
Vf = sum( 2 * pi * dr * ( r .* sf ) ); % numeric integration
if sf( end ) > sf( 1 ) % curving away
    Vf = pi * ( r( end )^2 - r(1)^2 ) * hf - Vf;
end
r = linspace( Db(1)/2, Db(2)/2, nr );
dr = ( Db(2) - Db(1) ) / 2 / nr;
sb = aspheric( [ -Rb(1) -Rb(1) kb ab ], r );
hb = sb( end );
Vb = sum( 2 * pi * dr * ( r .* sb ) ); % numeric integration
if sb( end ) > sb( 1 ) % curving away
    Vb = pi * ( r( end )^2 - r(1)^2 ) * hb - Vb;
end

Dmax = max( Df(2), Db(2) );
Dmin = min( Df(2), Db(2) );
R = Dmax / 2; % base radius
r = Dmin / 2; % top radius
Ve = pi/3 * h * ( R^2 + R * r + r^2 ); % conical volume at the equator
dmax = max( Df(1), Db(1) );
dmin = min( Df(1), Db(1) );
R = dmax / 2; % base radius for inside cut
r = dmin / 2; % top radius for inside cut
Ve = Ve - pi/3 * h * ( R^2 + R * r + r^2 ); % conical volume 0f the inside cut at the equator

ht = hf + hb + h; % total lens thickness in the center
V = Vf + Vb + Ve;

if disp_fl ~= 0 % create an animated gif
    % draw the lens
    bench = Bench();
    lensf = AsphericLens( [ -hf - h/2 0 0 ], Df, Rf, kf, af, { 'air', 'air' } );
    bench.append( lensf );
    % draw the outside wall
    if Df(2) == Db(2)
        lensh = CylinderLens( [ -h/2 0 0 ], Df(2), h, { 'air', 'air' } );
    else
        lensh = ConeLens( [ -h/2 0 0 ], Df(2), h, atan( ( Db(2) - Df(2) ) / 2 / h ), { 'air', 'air' } );
    end
    bench.append( lensh );
    if Df(1) ~= 0 || Db(1) ~= 0 % draw the inside wall
        hi = ht - sf(1) - sb(1);
        if Df(1) == Db(1)
            lensi = CylinderLens( [ -hf - h/2 + sf(1) 0 0 ], Df(1), hi, { 'air', 'air' } );
        else
            lensi = ConeLens( [ -hf - h/2 + sf(1) 0 0 ], Df(1), hi, atan( ( Db(1) - Df(1) ) / 2 / hi ), { 'air', 'air' } );
        end
        bench.append( lensi );
    end
    lensb = AsphericLens( [  h/2 + hb 0 0 ], Db, Rb, kb, ab, { 'air', 'air' } );
    bench.append( lensb );
    bench.draw( [], [], 0.33 ); % draw the lens with opaqueness set to 0.5
    camlight( 'left' );
    camlight( 'right' );
    camlight( 'headlight' );
    view( 0, 0 );
    
    filename = 'lens_dims.gif';
    for i = 1 : 36
        camorbit( 10, 0 );
        frame = getframe(1);
        im = frame2im( frame );
        [ imind, cm ] = rgb2ind( im, 256 );
        if i == 1
            imwrite( imind, cm, filename, 'gif', 'Loopcount', inf );
        else
            imwrite( imind, cm, filename, 'gif', 'WriteMode', 'append' );
        end
        pause( .02 );
    end
end

end

