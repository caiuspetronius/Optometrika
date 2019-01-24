classdef CylinderLens < Surface
    % CYLINDERLENS Implements a cylindrical lens surface.
    %
    % Member functions:
    %
    % l = CylinderLens( r, D, h, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - cylinder diameter, [ Dy Dz ] vector for an elliptical cylinder
    % h - cylinder height
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
        D = 1;    % cylinder diameter, [ Dy Dz ]
        h = 1;    % cylinder height
     end
    
    methods
        function self = CylinderLens( ar, aD, ah, aglass )
            if nargin == 0
                return;
            end
            self.D = aD;
            self.r = ar;
            self.h = ah;
            self.glass = aglass;
            self.R = aD / 2;  % this is needed for the correct handling of elliptic cylinders
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            if length( self.D ) == 1
                fprintf( 'Diameter:\t %.3f\n', self.D );
            else
                fprintf( 'Y Diameter:\t %.3f\n', self.D(1) );
                fprintf( 'Z Diameter:\t %.3f\n', self.D(2) );
            end
            fprintf( 'Height:\t %.3f\n', self.h );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self, color )
            % DISPLAY the cylindrical lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            nang = 100;
            [ x, y, z ] = cylinder( self.D(1) / 2, nang );
            z( 2, : ) = self.h;
            S = [ z(:) -y(:) x(:) ]; % put the cone into the Optometrika reference frame
            if length( self.D ) > 1  % elliptical cylinder
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
        end
        
    end
    
end

