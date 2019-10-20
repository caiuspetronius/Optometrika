function export_stl( lens1, lens2, centerThickness, nang, flange, filename )
% EXPORT_LENS exports the lens shape in STL format, requires stlwrite
% installed in the Matlab path
%
% INPUT:
%  lens1 - first lens surface lens object
%  lens2 - second lens surface lens object
%  centerThickness - thickness of the lens along the axis
%  nang - number of angular segments in a circle
%  flange - a 1 x 2 vector of flange heights
%  filename - name of the exported STL file
%
% Copyright: Yury Petrov, 2018
%

if nargin < 4
    nang = 100;
end
if nargin < 5
    flange = [];
end
if nargin < 6 || isempty( filename ) || strcmp( filename, '' )
    filename = [ datestr( now, 'mmddyyyyHHMMSS' ) '.stl' ];
end

as = linspace( 0, 2 * pi, nang );

%get lens profile
npoints = 1000;
if isprop( lens1, 'ncones' ) || isprop( lens2, 'ncones' ) % a Fresnel structure
    npoints = 5000;
    % npoints = round( D(1) / 2 * 1000 ); % to get radius steps in microns 
end

Y1 = linspace( lens1.D(1)/2, lens1.D(2)/2, npoints ); % starting from the center
Y2 = linspace( lens2.D(1)/2, lens2.D(2)/2, npoints ); % starting from the center

if isprop( lens1, 'R' ) && ~isempty( lens1.R )
    R1 = lens1.R;
    if length( R1 ) < 2
        R1 = [ R1 R1 ];
    end
    pars1 = [ R1(1) R1(2) lens1.k(1) ];
    if isprop( lens1, 'avec' ) && ~isempty( lens1.avec )
        pars1 = [ pars1 lens1.avec( :, 1 )' ];
    end
else
    pars1 = [ Inf Inf 0 ]; % planar surface 
end

if isprop( lens2, 'R' ) && ~isempty( lens2.R )
    R2 = lens2.R;
    if length( R2 ) < 2
        R2 = [ R2 R2 ];
    end
    pars2 = [ R2(1) R2(2) lens2.k(1) ];
    if isprop( lens2, 'avec' ) && ~isempty( lens2.avec )
        pars2 = [ pars2 lens2.avec( :, 1 )' ];
    end
else
    pars2 = [ Inf Inf 0 ]; % planar surface 
end

r = [ 0 centerThickness ];
for i = 1 : nang
    cs = cos( as( i ) );
    sn = sin( as( i ) );
    
    if ~isprop( lens1, 'ncones' ) % a non-Fresnel lens
        y1 = Y1;
        if ~isprop( lens1, 'funch' ) % not a general lens
            if isinf( pars1(1) ) && isinf( pars1(2) ) % a plane
                p1 = zeros( size( Y1 ) );
            else
                p1 = asphlens( Y1 * cs, Y1 * sn, pars1, 0 );
            end
        else
            p1 = lens1.funch( Y1 * cs, Y1 * sn, lens1.funca, 0 );
        end
    else % a Fresnel lens surface
        [ p1, y1 ] = frenlens( Y1, lens1 );
    end
    
    if ~isprop( lens2, 'ncones' ) % a non-Fresnel lens
        y2 = Y2;
        if ~isprop( lens2, 'funch' ) % not a general lens
            if isinf( pars2(1) ) && isinf( pars2(2) ) % a plane
                p2 = zeros( size( Y2 ) );
            else
                p2 = asphlens( Y2 * cs, Y2 * sn, pars2, 0 );
            end
        else
            p2 = lens1.funch( Y2 * cs, Y2 * sn, lens2.funca, 0 );
        end
    else % a Fresnel lens surface
        [ p2, y2 ] = frenlens( Y2, lens2 );
    end
    
    if isempty( flange )
        tmp = [ -fliplr( y1 ), y1, fliplr( y2 ), -y2, -y1( end ) ];
        z( i, : ) = tmp * sn;
        y( i, : ) = tmp * cs;
        x( i, : ) = [ r( 1 ) + fliplr( p1 ), r( 1 ) + p1, r( 2 ) + fliplr( p2 ), r( 2 ) + p2, r( 1 ) + p1( end ) ];
    else
        tmp = [ -fliplr( y1 ), y1, y1( end ) + flange(1), y2( end ) + flange(2), fliplr( y2 ), -y2, -y2( end ) - flange(2), -y1(end) - flange(1), -y1( end ) ];
        z( i, : ) = tmp * sn;
        y( i, : ) = tmp * cs;
        x( i, : ) = [ r(1) + fliplr( p1 ), r(1) + p1, r(1) + p1( end ), r(2) + p2( end ), r(2) + fliplr( p2 ), r(2) + p2, r(2) + p2( end ), r(1) + p1( end ), r(1) + p1( end ) ];
    end
end

%draw & export
figure;
h = surf( z, y, x, 'LineStyle','none', 'FaceColor', 'interp', 'FaceAlpha', 0.5 );
axis vis3d equal;

% export stl of the lens profile
if exist( 'stlwrite', 'file')
    [ f, v ] = surf2patch( h, 'triangles' );
    TP = triangulation( f, v );
    stlwrite( TP, filename );  % this is a build-in function in Matlab 2018b
    fprintf( 'Exported to %s\n', filename );
end
