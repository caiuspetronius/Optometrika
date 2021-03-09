classdef FresnelLens < Surface
    % FRESNELLENS Implements a Fresnel lens surface given by a rotation of
    % a piecewize linear curve defining a set of cone segments. These are 
    % characterized by radii, sag points, and half-angles measured with
    % respect to the lens axis of rotation.
    %
    % Member functions:
    %
    % l = FresnelLens( r, rad, sag, the, glass, walls, inner ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % rad - ncones x 1 vector of outer radii for the Fresnel cones
    % sag - ncones x 1 vector of sag values for the inner cone edge
    % the - ncones x 3 matrix of cone half-angles (radians) and curvature
    % radius (in units of the outer ring radius).
    % If the second column is absent, its half angles are assuned
    % to be equal to pi, i.e., the Fresnel walls - vertical cylinders.
    % If the third column is absent, it's assumed to be Inf, i.e., flat
    % conical rings.
    % Values in the second column less than pi indicate conical walls sloping
    % toward the lens periphery, the values more than pi - conical walls
    % sloping toward the lens center (overhanging).
    % glass - 1 x 2 cell array of strings, e.g., { 'air' 'acrylic' }
    % walls - 0 for soot Fresnel walls (absorbs all light), 1 for walls of the same
    % material as the rest of the lens
    % inner - 0 for inner Fresnel diameter = 0, 1 for inner Fresnel
    % diameter = 2 * ( 2 * rad(1) - rad(2) ), i.e. defined by the smallest
    % inner radius of the Fresnel ring
    %
    % OUTPUT:
    % l - lens surface object
    %
    % l.display() - displays the surface l information
    %
    % l.draw() - draws the surface l in the current axes
    %
    % l.rotate( rot_axis, rot_angle ) - rotate the surface l
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    % 
    % Copyright: Yury Petrov, 2017
    %
    
    properties
        D = [ 0; 1 ]; % lens diameter (inner, outer)
        ncones = 0; % number of Fresnel cones
        rad = [];   % cone radii
        sag = [];   % cone sags
        the = [];  % cone half-angles
        vx = [];   % conic vertex (tip of ellipsoid, parabola, hyperbola) x position
        walls = 1; % consider Fresnel walls by default
        possible = true; % if the lens is possible or not
     end
    
    methods
        function self = FresnelLens( ar, aR, as, ath, aglass, awalls, inner )
            if nargin == 0
                return;
            end
            if nargin < 6
                awalls = 1; % consider Fresnel walls by default
            end
            
            if size( aR, 1 ) < size( aR, 2 )
                aR = aR';
            end
            if size( as, 1 ) < size( as, 2 )
                as = as';
            end
            if size( ath, 1 ) < size( ath, 2 )
                ath = ath';
            end
            if size( aR, 1 ) ~= size( as, 1 ) || ...
               size( aR, 1 ) ~= size( ath, 1 ) || ...
               size( as, 1 ) ~= size( ath, 1 )
                error( 'Fresnel lens R, s, and th vectors must have the same length!' );
            else
                self.ncones = size( aR, 1 );
            end
            % sort all vectors based on ascending radii
            [ aR, ind ] = sort( aR );
            as = as( ind );
            ath = ath( ind, : );
            if size( ath, 2 ) == 1
                ath = [ ath pi * ones( size( ath ) ) ]; % default to vertical walls
            end
            
            if nargin < 7 || inner == 0
                self.D(1) = 0;
            else
                self.D(1) = 2 * ( aR( 1 ) - ( aR(2) - aR(1) ) ); % inner radius, assuming first and second rings have the same width
            end
            self.D(2) = 2 * aR( end ); % outer radius
           
            self.r = ar;
            self.sag = as;
            self.the = ath( :, 1:2 );
            self.walls = awalls;
            
            if size( ath, 2 ) == 3 % curvature specified
                self.R = ath( :, 3 ) .* aR; % radius of curvature
                rs = [ self.D(1)/2; aR ];
                rs = rs( 1 : end - 1 ); % inner radii
                self.k = ( self.R ./ rs ).^2 - tan( self.the( :, 1 ) ).^2 - 1; % conic constants
                self.vx = -self.R ./ ( 1 + self.k ) .* ( 1 - tan( self.the( :, 1 ) ) .* rs ./ self.R ); % quadric vertices using the inner radius
                if rs( 1 ) == 0  % k cannot be defined by the slope at the center of the quadric, take it explicitely
                    self.k( 1 ) = self.the( 1, 1 );
                    self.the( 1 ) = pi/2;
                    self.vx( 1 ) = 0; % because the inner radius is zero here
                end
                % treat parabolas specially
                parabola = abs( self.k + 1 ) < 1e-10;
                ellipse = self.k + 1 >= 1e-10;
                self.k( parabola ) = -1; % round to avoid errors at the ray-tracing step
                self.vx( parabola ) = -rs( parabola ).^2 ./ ( 2 * self.R( parabola ) );
                if parabola( 1 ) && rs( 1 ) == 0
                    self.vx( 1 ) = -rs( 2 ).^2 ./ ( 2 * self.R( 1 ) );
                end
                    
                % find intersection of the quadric with the wall cone
                sl = tan( pi/2 - self.the( :, 2 ) ); % slope of the wall
                g = 1 ./ ( 1 + self.k );
                dd = circshift( self.sag, -1 ) - self.sag;
                dd( end ) = dd( end - 1 ); % approximate with linear here
                z0 = dd - self.vx;
                a = sl.^2 + g;
                b = 2 * sl .* ( z0 - self.R .* g - aR .* sl );
                c = sl.^2 .* aR.^2 - 2 * sl .* aR .* ( z0 - self.R .* g ) + z0.^2 - 2 * z0 .* self.R .* g;
                % treat parabolas specially
                a( parabola ) = 1 ./ ( 2 * self.R( parabola ) );
                b( parabola ) = -sl( parabola );
                c( parabola ) = sl( parabola ) .* aR( parabola ) - z0( parabola );
                % solve quadratic equation
                D2 = b.^2 - 4 * a .* c;
                D = sqrt( D2 );
                sgn = 2 * ( self.the( :, 2 ) > pi ) - 1;
                sgn( parabola | ellipse ) = -sgn( parabola | ellipse );
                rm = 0.5 * ( -b + sgn .* sign( self.R ) .* D ) ./ a; % use the smaller of the two intersections
                dd = rm - circshift( aR, 1 ); % distance between the wall intersection and the previous outer radius
                if sum( dd( 2 : end ) < 0 ) ~= 0 % the quadric curves back
                    self.possible = false;
                end
                % rm( end ) = aR( end ); % vertical wall for the last ring
                aR = [ aR rm ];      
                % deal with vertical walls (cylinders)
                vertical = abs( self.the( :, 2 ) - pi ) < 1e-10; 
                aR( vertical, 2 ) = aR( vertical, 1 ); % the same radius here
                % check if the intersections are possible
                if sum( D2( ~vertical ) < 0 ) > 0
                    self.possible = false;
                end
                           
            else % conic rings
                % find intersection of the Fresnel cone with the wall cone
                if size( ath, 2 ) > 1 % non-cylindrical (conic) walls
                    %tga = tan( pi - ath( :, 1 ) );
                    %tth = tan( ath( :, 2 ) );
                    tga = tan( ath( :, 1 ) );
                    tth = tan( pi - ath( :, 2 ) );
                    dr = aR - circshift( aR, 1 );
                    dr( 1 ) = aR( 1 );
                    ds = as - circshift( as, 1 );
                    ds( 1 ) = as( 1 );
                    dx = -tth .* ( dr ./ tga - ds ) ./ ( tth ./ tga + 1 ); % distance from the next radius aR to the intersection of the Fresnel slope and wall
                    aR = [ aR, aR + dx ]; % the first radius is the value where the next Fresnel slope starts,
                    % the second value is where the current Fresnel slope ends and the Fresnel wall starts
                else
                    aR = [ aR, aR ];
                    ath = [ ath, repmat( pi, size( ath, 1 ), 1 ) ];
                end               
            end
            
            self.rad = aR;
            self.glass = aglass;
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D(2) );
            if self.D(1) ~= 0
                fprintf( 'Inner diameter:\t %.3f\n', self.D(1) );
            end
            fprintf( 'Number of Fresnel rings:\t %i\n', self.ncones );
            fprintf( 'Steepest slope (rad):\t %i\n', max( abs( self.the( :, 1 ) ) ) );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self, color )
            % DISPLAY the Fresnel lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            nang = 100;
            as = linspace( 0, 2 * pi, nang );
            nrad = 100;
            h = zeros( self.ncones, 1 );
            cnt = 0; % surface counter
            
%             F = [];
%             V = [];
            for i = 1 : self.ncones % loop over cones
                zc = [];
                if i == 1
                    if self.D(1) == 0
                        radin = 0; 
                    else
                        radin = 2 * self.rad( 1, 1 ) - self.rad( 2, 1 ); % take the inner radius assuming the same step
                    end
                else
                    radin = self.rad( i - 1, 1 );
                end
                
                if isempty( self.vx ) || self.R( i, 1 ) == 0 || isinf( self.k( i ) ) % cone surface
                    [ x, y, z ] = cylinder( [ radin self.rad( i, 2 ) ], nang );
                    z( 1, : ) = self.sag( i );
                    z( 2, : ) = self.sag( i ) + ( self.rad( i, 2 ) - radin ) / tan( self.the( i, 1 ) );
                    [ xc, yc ]  = cylinder( [ self.rad( i, 2 ) self.rad( i, 1 ) ], nang ); % Fresnel wall
                else % quadric surface                  
                    rads = linspace( radin, self.rad( i, 2 ), nrad );
                    [ angs, rads ] = meshgrid( as, rads );
                    [ x, y ] = pol2cart( angs, rads );
                    a = 1 + self.k( i );
                    if length( self.R( i, : ) ) == 1
                        r2xy = ( x.^2 + y.^2 ) / self.R( i )^2;
                        if abs( a ) < 1e-10 % paraboloid, special case
                            z = r2xy * self.R( i, 1 ) / 2;
                        else
                            z = self.R( i ) / a * ( 1 - sqrt( 1 - a * r2xy ) );
%                             if sum( ~isreal( z ) ) ~= 0 
%                                 a;
%                             end
                        end
                    else % asymmetric conic
                        r2xy = x.^2 / self.R( i, 1 ) + y.^2 / self.R( i, 2 );
                        if abs( a ) < 1e-10 % paraboloid, special case
                            z = r2xy / 2;
                        else
                            z = r2xy ./ ( 1 + sqrt( 1 - a * ( x.^2 / self.R( i, 1 )^2 + self.R( i, 2 ) / self.R( i, 1 ) * y.^2 / self.R( i, 2 )^2 ) ) );
                        end
                    end
                    z = z + self.sag( i ) + self.vx( i ); % add sag and quadric vertex displacement
                    [ xc, yc ]  = cylinder( [ self.rad( i, 2 ) self.rad( i, 1 ) ], nang - 1 ); % Fresnel wall
                end
                
                % create the cylinder wall
                zc( 1, : ) = z( end, : );
                if i == self.ncones
                    zc( 2, : ) = zc( 1, : ); % no wall for the last ring
                else
                    zc( 2, : ) = self.sag( i + 1 ); % add sag to the first point
                end
                
                if abs( self.the( i, 1 ) - pi/2 ) > 1e-10 && ( sign( ( z( 2, 1 ) - z( 1, 1 ) ) / ( self.rad( i, 2 ) - radin ) ) ~= sign( tan( self.the( i, 1 ) ) ) )
                    z = flipud( z );
                end
                if sign( ( zc( 2, 1 ) - zc( 1, 1 ) ) / ( self.rad( i, 1 ) - self.rad( i, 2 ) ) ) ~= sign( tan( self.the( i, 2 ) ) )
                    zc = flipud( zc );
                end
                S  = [ z(:) -y(:)   x(:) ]; % put the cone into the Optometrika reference frame
                Sc = [ zc(:) -yc(:) xc(:) ]; % put the wall into the Optometrika reference frame
                
                % rotate and shift
                if self.rotang ~= 0
                    S =  rodrigues_rot( S,  self.rotax, self.rotang );
                    Sc = rodrigues_rot( Sc, self.rotax, self.rotang );
                end
                x(:) =  S( :, 1 ) + self.r( 1 );
                y(:) =  S( :, 2 ) + self.r( 2 );
                z(:) =  S( :, 3 ) + self.r( 3 );
                xc(:) = Sc( :, 1 ) + self.r( 1 );
                yc(:) = Sc( :, 2 ) + self.r( 2 );
                zc(:) =  Sc( :, 3 ) + self.r( 3 );
                
                c = repmat( reshape( color( 1:3 ), [ 1 1 3 ] ), size( x, 1 ), size( x, 2 ), 1 );
                cnt = cnt + 1;
                if sum( ~isreal( x ) + ~isreal( y ) + ~isreal( z ) ) ~= 0
                    error( 'Complex surface!' );
                end
                h( cnt ) = surf( x, y, z, c, ...
                          'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', color(4), ...
                          'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
                fv = surf2patch( x, y, z, 'triangles' );
                %F = [ F; fv.faces ];
                %V = [ V; fv.vertices ];
                c = repmat( reshape( color( 1:3 ), [ 1 1 3 ] ), size( xc, 1 ), size( xc, 2 ), 1 );
                cnt = cnt + 1;
                h( cnt ) = surf( xc, yc, zc, c, ...
                          'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', color(4), ...
                          'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
                fv = surf2patch( xc, yc, zc, 'triangles' );
                %F = [ F; fv.faces ];
                %V = [ V; fv.vertices ];
            end
             
            %stlwrite( 'test.stl', F, V );
        end 
        
        function rotate( self, rot_axis, rot_angle )
            self.rotate@Surface( rot_axis, rot_angle ); % rotate the surface members
            if abs( rot_angle ) > pi/2
                self.the = pi - self.the;
                self.the( self.the == 0 ) = pi;
                self.sag = -self.sag;
                self.vx  = -self.vx;
            end
        end
        
    end
    
end

