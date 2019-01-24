classdef ConeLens < Surface
    % CONELENS Implements a cone lens surface.
    %
    % Member functions:
    %
    % l = ConeLens( r, D, the, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - 2 x 1 (inner outer diameters) or 1 x 1 vector (outer diameter)
    % the - cone half-angle, radians
    % glass - 1 x 2 cell array of strings, e.g., { 'air' 'acrylic' }
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
    % Copyright: Yury Petrov, 2016
    %
    
    properties
        D = [ 0; 1 ];     % lens diameter: inner/outer
        rad = [ 0; 1 ]   % cone radii
        the = pi/4;      % cone half-angle
     end
    
    methods
        function self = ConeLens( ar, aD, ath, aglass )
            if nargin == 0
                return;
            end
            if size( aD, 1 ) < size( aD, 2 )
                aD = aD';
            end
            if size( aD, 1 ) == 1
                aD = [ 0; aD ]; % assume inner radius = 0
            end        
            self.D = aD;
            self.r = ar;
            self.rad = aD / 2;
            self.the = ath;
            self.glass = aglass;
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D(2) );
            if self.D(1) ~= 0
                fprintf( 'Inner diameter:\t %.3f\n', self.D(1) );
            end
            fprintf( 'Slope (rad):\t %.3f\n', abs( self.the ) );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self, color )
            % DISPLAY the cone lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            nang = 100;
            [ x, y, z ] = cylinder( self.rad, nang );
            z( 2, : ) = ( self.rad( 2 ) - self.rad( 1 ) ) / tan( self.the );
            S = [ z(:) -y(:) x(:) ]; % put the cone into the Optometrika reference frame
            
            % rotate and shift
            if self.rotang ~= 0
                S = rodrigues_rot( S, self.rotax, self.rotang );
            end
            x(:) = S( :, 1 ) + self.r( 1 );
            y(:) = S( :, 2 ) + self.r( 2 );
            z(:) = S( :, 3 ) + self.r( 3 );
            
            c = repmat( reshape( color( 1:3 ), [ 1 1 3 ] ), size( x, 1 ), size( x, 2 ), 1 );
            h = surf( x, y, z, c, ...
                'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', color(4), ...
                'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
        end
        
        function rotate( self, rot_axis, rot_angle )
            self.rotate@Surface( rot_axis, rot_angle ); % rotate the surface members
            if abs( rot_angle ) > pi/2
                self.the = pi - self.the;
            end
        end
        
    end
    
end

