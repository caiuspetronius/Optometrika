classdef ConeLens < Surface
    % CONELENS Implements a cone lens surface.
    %
    % Member functions:
    %
    % l = ConeLens( r, D, h, the, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - cone top (smaller x coordinate) diameter, [ Dy Dz ] vector for an elliptical cone
    % h - cone height
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
        D = 0;           % cone top diameter
        h = 1;           % cone height
        rad = [ 0; 1 ]   % cone top/bottom radii
        the = pi/4;      % cone half-angle
     end
    
    methods
        function self = ConeLens( ar, aD, ah, ath, aglass )
            if nargin == 0
                return;
            end
            self.r = ar;
            self.D = aD;
            self.h = ah;
            self.the = ath;
            self.glass = aglass;
            self.rad(1) = aD(1) / 2; % top radius along y-axis
            self.rad(2) = self.rad(1) + ah * tan( ath ); % bottom radius along y-axis
            self.R = aD / 2;    % this is needed for the correct handling of elliptic cones
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            if length( self.D ) == 1
            fprintf( 'Diameter:\t %.3f\n', self.D );
            else % elliptic cylinder
                fprintf( 'Y Diameter:\t %.3f\n', self.D(1) );
                fprintf( 'Z Diameter:\t %.3f\n', self.D(2) );
            end
            fprintf( 'Height:\t %.3f\n', self.h );
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
            z( 2, : ) = self.h;
            S = [ z(:) -y(:) x(:) ]; % put the cone into the Optometrika reference frame
            if length( self.D ) > 1  % elliptical cone
                S( :, 3 ) = S( :, 3 ) * self.D(2) / self.D(1);
            end
            
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

