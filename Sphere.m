 classdef Sphere < Surface
    % SPHERE Implements a surface as a sphere

    %
    % Member functions:
    %
    % l = Lens( r, D, R, k, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - diameter, 1x1 vector (outer)
    % k - conic coefficient, for astigmatic surface corresponds to the y-axis 
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
        D = [ 1 ]; % lens diameter (inner, outer)
     end
    
    methods
        function self = Sphere(ar, aD, aglass)
            if nargin == 0
                return;
            end
            self.r = ar;
            self.D = aD;
            self.glass = aglass;
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D(2) );
            if self.D(1) ~= 0
                fprintf( 'Inner diameter:\t %.3f\n', self.D(1) );
            end
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self, color )
            % DISPLAY the lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end      
            [X,Y,Z] = sphere(100);
            clear x y z;
            x(:,:) = X(:,:)*self.D/2 + self.r( 1 );
            y(:,:) = Y(:,:)*self.D/2 + self.r( 2 );
            z(:,:) = Z(:,:)*self.D/2 + self.r( 3 );
            
            c = repmat( reshape( color( 1:3 ), [ 1 1 3 ] ), size( x, 1 ), size( x, 2 ), 1 );
            h = surf( x, y, z, c, ...
                'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', color(4), ...
                'AmbientStrength', 0., 'SpecularStrength', 1 ); % grey color, shiny
        end
        
    end
    
end

