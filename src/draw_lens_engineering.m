function [ f, pr ] = draw_lens_engineering( r, D, R, k, a, glass, name, varargin  )
%
% Make a pdf figure with the engineering drawing of a lens. Also exports
% DXF file of its profile, if dxflib are installed in the MATLAB path.
%
% INPUT:
%   r - 1 x 2 vector of surface locations along the x-axis (optical axis)
%   D - 1 x 2 vector of the lens surface diameters
%   R - 1 x 2 vector of the lens surface radii of curvature
%   k - 1 x 2 vector of the lens conic constants
%   a - n x 2 vector of aspherical polynomial terms, one column per surface
%   glass - string with the lens material name
%   name - name of the drawing displayed as its title, defaults to current
%          date and time
%   varargin - strings specifying other manufacturing options. 'AR' adds
%          instructions to use broadband antireflective coating, 'flange' 
%          followed by 2 numbers specifies a flange heights, other options 
%          are up to you to specify.
%
% OUTPUT:
%   f - figure handle
%   pr - a cell array with high-precision profiles of the two lens surfaces
%
% EXAMPLE:
%   f = draw_lens_engineering( [ 0 12 ], [ 40 40 ], [ 50 -30 ], [ -2 -1 ], [], 'zeonex', [], 'AR', 'Make my lens really beautiful')
%
% Copyright: Yury Petrov, 2017
%


if nargin < 7 || isempty( name ) || strcmp( name, '' )
    name = datestr( now, 'mmddyyyyHHMMSS' );
end
if nargin < 6
    error( 'Please specify at least the first 6 arguments!' );
end
if size( a, 2 ) > 2
    a = a';
    if size( a, 2 ) ~= 2
        error( 'Matrix of polynomial coefficients has to be n x 2!' );
    end
end
if sum( abs( a(:) ) ) == 0
    a = [];
end

% get screen size 
ScrSz = get( 0, 'ScreenSize' ); % desired aspect ratio 3:4 -> check and adjust target figure size 
if ScrSz(3) / ( 0.95 * ScrSz(4) ) > 4/3 %vertical space limits (FHD etc.) 
    figh = 0.95 * ScrSz(4); %leave some space for taskbars, which could have varying sizes in different OS 
    figw = figh * 4/3; 
else %horizontal space limits (could be in portrait mode) 
    figw = ScrSz(3); 
    figh = figw * 3/4; 
end
f = figure( 'Name', name, 'Position', [ 0 0 figw figh ], 'Color', [ 1 1 1 ] );

% title
an = annotation( 'textbox', [ 0.65 0.02 0.33 0.08 ], 'String', 'Title', 'FitBoxToText', 'off' );
an.FontSize = 18;
an.Units = 'normalized';
an = text( 0.75, -0.06, name, 'FontSize', 36 ); 
an.Units = 'normalized';
axis off;

% table with lens parameters
t = uitable( f );

t.ColumnName = { '' 'S1' 'S2' };
t.ColumnWidth = { 200 120 120 };
t.ColumnFormat = { 'char' 'char' 'char' };
t.RowName = {};
t.RowStriping = 'off';
t.FontSize = 18;

s1 = 'Convex';
if R(1) < 0
    s1 = 'Concave';
end
s2 = 'Convex';
if R(2) > 0
    s2 = 'Concave';
end

t.Data = { 
    'Shape', s1, s2;
    'Radius', R(1), R(2);
    'Aspherical conic', k(1), k(2);
    };
if ~isempty( a ) && ~sum( abs( a(:) ) ) == 0
    for i = 1 : size( a, 1 )
        t.Data = [ t.Data; { [ 'Polynomial r^' num2str( i * 2 ) ], a( i, 1 ), a( i, 2 ); } ];
    end
end
t.Data = [ t.Data; 
      {'Clear aperture', D(1), D(2);
       'Surface figure error', '6 fringes', '6 fringes';
       'Surface irregularity', '2 fringes', '2 fringes';
       'Surface roughness', '10 nm', '10 nm';
       'Scratch-dig', '40-20', '40-20' 
      }
    ];

margin = 40; % margin for the table position, px
parentPosition = get( f, 'position' );
t.Position( 3:4 ) = t.Extent( 3:4 );
%t.Position( 1:2 ) = [ parentPosition(1) + margin, parentPosition(2) + parentPosition(4) - t.Position(4) - margin ];
t.Position( 1:2 ) = [ parentPosition(1) + margin, parentPosition(4) - t.Position(4) - margin ];

% text with other specs
Fr = {}; % empty cell to store Fresnel structures
Frcnt = 0; % Fresnel surfaces counter
str = { 
    '1. All dimensions in millimeters', ...
    [ '2. Allowable optical material: ' glass], ... 
    '3. Center element within 0.005 mm TIR', ...
    '4. Lens wedge: 0.02 mm ETD', ...
    '5. Fine ground on lens edge surface', ...
    '6. Lens fabrication via diamond turning' };
cnt = 6;
flange = [];
for i = 1 : length( varargin )
    if isempty( varargin{ i } )
        continue;
    end
    cnt = cnt + 1;
    if strcmpi( varargin{ i }, 'AR' )
        str = [ str, ...
            [ num2str( cnt ) '. BBAR coating on both surfaces' ], ...
            '    Ravg < 1%', ...
            '    Spectrum range: 450 - 650 nm', ...
            '    Angle of incidence: 0 - 45 degrees' ];
    elseif strcmpi( varargin{ i }, 'flange' )
        flange = [ varargin{ i + 1 } varargin{ i + 2 } ];
        if ~isnumeric( flange )
            error( 'Two flange heights must be numeric values!' );
        end
        varargin{ i + 1 } = []; % skip the next two values
        varargin{ i + 2 } = []; % skip the next two values
        str = [ str, [ num2str( cnt ) '. Flange height [' num2str( flange(1) ) ', ' num2str( flange(2) ) ']' ] ];
    elseif strcmpi( varargin{ i }, 'Fresnel' )
        Frcnt = Frcnt + 1;
        Fr{ Frcnt } = varargin{ i + 1 }; % read the Fresnel structure
        varargin{ i + 1 } = []; % skip the next value
    else
        str = [ str, [ num2str( cnt ) '. ' varargin{ i } ] ];
    end
end
m = uicontrol( f, 'Style', 'text', 'String', str, 'Position', [ margin t.Position(2) - 450 400 400 ], 'BackgroundColor', [ 1 1 1 ], ...
    'FontSize', 18, 'HorizontalAlignment', 'left' );

% draw the lens
clr = [ 0 .3 .6 ]; % dark blue
if isempty( Fr ) % both surfaces regular
    lw = 2; % line width
else
    lw = 1;
end

% draw the lens front view (circles)
a1 = axes( 'Position', [ 0.4 0.3 0.3 0.3 ], 'Clipping', 'off' );
if ~isempty( flange )
    circle( [ D, D + flange ], [ r(1), r(2), 1.001 * r(1), 0.999 * r(2) ], clr, lw ); % make flange edge so slightly behind the front face edge and
else
    circle( D, [ r(1), r(2) ], clr, lw );
end
% in front of the back face edge to ensure right depth masking
axis equal tight off;

% draw the lens profile view
a2 = axes( 'Position', [ 0.7 0.3 0.2 0.3 ], 'Clipping', 'off' );
npoints = 1000;
if ~isempty( Fr ) > 0 && ~isempty( Fr{1} ) % a Fresnel structure
    npoints = 10000;
    % npoints = round( D(1) / 2 * 1000 ); % to get radius steps in microns 
end
z1 = zeros( 1, npoints );
y1 = linspace( 0, D(1)/2, npoints ); % starting from the center
npoints = 1000;
if length( Fr ) > 1 && ~isempty( Fr{2} ) % a Fresnel structure
    npoints = 10000;
    % npoints = round( D(2) / 2 * 1000 ); % to get radius steps in microns 
end
z2 = zeros( 1, npoints );
y2 = linspace( 0, D( 2 )/2, npoints ); % starting from the center
if isempty( a )
    if isempty( Fr ) % a regular conic lens
        p1 = asphlens( y1, z1, [ R( 1 ) R( 1 ) k( 1 ) ], 0 );
        p2 = asphlens( y2, z2, [ R( 2 ) R( 2 ) k( 2 ) ], 0 );
    else % a Fresnel lens
        if length( Fr ) == 1
            [ p1, y1 ] = frenlens( y1, Fr{1} );
            p2 = asphlens( y2, z2, [ R( 2 ) R( 2 ) k( 2 ) ], 0 );
        else
            if isempty( Fr{ 1 } ) % back surface Fresnel
                p1 = asphlens( y1, z1, [ R( 1 ) R( 1 ) k( 1 ) ], 0 );
                [ p2, y2 ] = frenlens( y2, Fr{2} );
            else % both surfaces Fresnel
                [ p1, y1 ] = frenlens( y1, Fr{1} );
                [ p2, y2 ] = frenlens( y2, Fr{2} );
            end
        end
    end 
else
    p1 = asphlens( y1, z1, [ R( 1 ) R( 1 ) k( 1 ) a( :, 1 )' ], 0 );
    p2 = asphlens( y2, z2, [ R( 2 ) R( 2 ) k( 2 ) a( :, 2 )' ], 0 );
end
if isempty( flange )
    y = [ -fliplr( y1 ), y1, fliplr( y2 ), -y2, -y1( end ) ];
    x = [ r( 1 ) + fliplr( p1 ), r( 1 ) + p1, r( 2 ) + fliplr( p2 ), r( 2 ) + p2, r( 1 ) + p1( end ) ];
else
    y = [ -fliplr( y1 ), y1, y1( end ) + flange(1), y2( end ) + flange(2), fliplr( y2 ), -y2, -y2( end ) - flange(2), -y1(end) - flange(1), -y1( end ) ];
    x = [ r(1) + fliplr( p1 ), r(1) + p1, r(1) + p1( end ), r(2) + p2( end ), r(2) + fliplr( p2 ), r(2) + p2, r(2) + p2( end ), r(1) + p1( end ), r(1) + p1( end ) ];
end

pr = { [ y1' p1' ] [ y2' p2' ] }; % radii and corresponding surface profiles

p = patch( x, y, [ 1 0 0 ], 'EdgeColor', clr, 'LineWidth', lw );
hatchfill( p, 'single', 30, 35, 'none', clr );
axis equal tight off;

if exist( 'dxf_open' ) % export dxf of the lens profile
    FID = dxf_open( [ name '_profile.dxf' ] );
    dxf_polyline( FID, x', y', zeros( size( x ) )' );
    dxf_close(FID);
end

% rads = y;
% z = x;
% nang = 100;
% as = linspace( 0, 2 * pi, nang );
% [ angs, rads ] = meshgrid( as, rads );
% [ x, y ] = pol2cart( angs, rads );
% figure;
% h = surf( x, y, repmat( z, nang, 1 )', 'LineStyle','none', 'FaceColor', 'interp', 'FaceAlpha', 0.5 ); axis vis3d equal;
% if exist( 'stlwrite' ) % export dxf of the lens profile
%     fv = surf2patch( h, 'triangles' );
%     stlwrite( [ name '.stl' ], fv );
% end

% draw the dimensions
R = max( D ) / 2;
Rm = min( D ) / 2;
an = annotation( 'line' );
set( an, 'Parent', a1 );
set( an, 'X', [ -R R ], 'Y', [ 0, 0 ], 'LineStyle', '-.', 'LineWidth', 1.5 ); % center horizontal diameter line
an = annotation( 'line' );
set( an, 'Parent', a1 );
set( an, 'X', [ 0 0 ], 'Y', [ -1.5 * Rm, 1.5 * Rm ], 'LineStyle', '-.', 'LineWidth', 1.5 ); % center vertical diameter line
for i = 0 : 1
    sgn = 1 - 2 * i;
    an = annotation( 'arrow' );
    set( an, 'Parent', a1 );
    set( an, 'X', [ 0, 0.4 * R ], 'Y', sgn * [ 1.5 * Rm, 1.5 * Rm ], 'LineWidth', 2 ); % projection arrow
    text( a1, -0.2 * R, sgn * 1.5 * Rm, 'A', 'FontSize', 30 );
end
ei = 2;
if D(1) == D(2)
    ei = 1;
end
for i = 1 : ei % loop over radii
    k = 1.75;
    if i == 2
        k = -k;
    end
    R = D(i)/2;
    an = annotation( 'line' );
    set( an, 'Parent', a1 );
    set( an, 'X', [ -R -R ], 'Y', [ 0, k * R ] ); % left diameter line
    an = annotation( 'line' );
    set( an, 'Parent', a1 );
    set( an, 'X', [  R  R ], 'Y', [ 0, k * R ] ); % right diameter line
    an = annotation( 'textarrow', 'String', [ '\oslash ' num2str( D( i ) ) ], 'FontSize', 18 );
    set( an, 'Parent', a1 );
    set( an, 'X', [  0.15 * R,  R ], 'Y', [ k * R, k * R ] ); % diameter dimension double arrow
    an = annotation( 'arrow' );
    set( an, 'Parent', a1 );
    set( an, 'X', [  -0.2 * R,  -R ], 'Y', [ k * R, k * R ] ); % diameter dimension double arrow
end

% cross-section dimensions
Rm = max( D ) / 2;
text( a2, mean( r ) - 1.2 * max( abs( [ p1 p2 ] ) ), 1.7 * Rm, 'Section A-A', 'FontSize', 30 );
k = 1.25; % overall dimensions margin scaling factor
al = 0.25; % arrow length scaling factor
ts = 0.1; % relative text shift wrt dimensions arrow
for i = 1 : 2
    x = sort( r );
    r0 = [ 0 0 ];
    if i == 2
        k = -k;
        x = [ r(1) + p1( end ), r(2) + p2( end ) ];
        r0 = -D / 2;
    end
    R = D(i)/2;
    an = annotation( 'line' );
    set( an, 'Parent', a2 );
    set( an, 'X', [ x(1) x(1) ], 'Y', [ r0(1), k * Rm ] ); % left thickness line
    an = annotation( 'line' );
    set( an, 'Parent', a2 );
    set( an, 'X', [ x(2) x(2) ], 'Y', [ r0(2), k * Rm ] ); % right thickness line
    if x(2) - x(1) > 0.5 * Rm % wide enough for double arrow
        an = annotation( 'doublearrow' );
        set( an, 'Parent', a2 );
        set( an, 'X', x, 'Y', [ k * Rm, k * Rm ] ); % diameter dimension double arrow
    else
        an = annotation( 'arrow' );
        set( an, 'Parent', a2 );
        set( an, 'X', [ x(1) - al * Rm, x(1) ], 'Y', [ k * Rm, k * Rm ] ); % diameter dimension double arrow
        an = annotation( 'arrow' );
        set( an, 'Parent', a2 );
        set( an, 'X', [ x(2) + al * Rm, x(2) ], 'Y', [ k * Rm, k * Rm ] ); % diameter dimension double arrow
    end
    text( a2, -0.1 * R + mean( x ), k * ( 1 + ts ) * Rm, num2str( x(2) - x(1), '%.2f' ), 'FontSize', 18 );
    
    % display flange dimensions
    g = 1.5 * al * ( 2 * i - 3 );
    if ~isempty( flange ) && ~( flange(1) == flange(2) && i == 1 ) % draw only one flange dimension if the flange is symmetric
        an = annotation( 'line' );
        set( an, 'Parent', a2 );
        x0 = r( i ) + pr{ i }( end );
        set( an, 'X', [ x0 + g * Rm, x0 ], 'Y', -[ R, R ] ); % top thickness line
        an = annotation( 'line' );
        set( an, 'Parent', a2 );
        set( an, 'X', [ x0 + g * Rm, x0 ], 'Y', -flange( i ) - [ R, R ] ); % bottom thickness line
        an = annotation( 'arrow' );
        set( an, 'Parent', a2 );
        set( an, 'X', [ x0 + g * Rm, x0 + g * Rm ], 'Y', -[ R - al * Rm, R ] ); % diameter dimension double arrow
        an = annotation( 'arrow' );
        set( an, 'Parent', a2 );
        set( an, 'X', [ x0 + g * Rm, x0 + g * Rm ], 'Y', -flange( i ) - [ R + al * Rm, R ] ); % diameter dimension double arrow
        text( a2, x0 + g * Rm + ( i == 2 ) * 0.5 * ts * Rm - ( i == 1 ) * 3 * ts * Rm, -R - flange( i )/2, num2str( flange( i ), '%3.2f' ), 'FontSize', 18 );
    end
end

% mark surfaces
text( a2, r(1) - 0.25 * Rm, 0, 'S1', 'FontSize', 24 );
text( a2, r(2) + 0.05 * Rm, 0, 'S2', 'FontSize', 24 );

% print figure to pdf setting the right size and resolution for the print
set( f, 'Units', 'inches' );
pos = get( f, 'Position');
set( f, 'PaperOrientation', 'landscape', 'PaperSize', [ pos(3) pos(4) ] );
print( f, [ name '.pdf' ], '-dpdf', '-r0' );



function circle( D, x, color, lw )
% draw concentric circle(s) of a given diameter and color at x depths
[ ~, is ] = sort( x ); % sort by depth (increasing)
D = D( is ); % order diameters
for i = 1 : length( D )
    pos = [ -D( i )/2,  -D( i )/2,  D( i ), D( i ) ];
    if sum( D( 1 : i - 1 ) > D( i ) ) == 0 % all previous diameters were smaller
        rectangle( 'Position', pos, 'Curvature', [ 1 1 ], 'EdgeColor', color, 'LineWidth', lw ); % draw solid circle
    else
        rectangle( 'Position', pos, 'Curvature', [ 1 1 ], 'EdgeColor', color, 'LineWidth', 1, 'LineStyle', '--' ); % draw dashed circle
    end
end
axis equal;
