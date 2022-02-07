 classdef Lens < Surface
    % LENS Implements a lens surface given by a rotation of a conic curve
    % (conic) lens surface given by
    % z = 1/R * r^2 ./ ( 1 + sqrt( 1 - ( 1 + k ) * (r/R)^2 ) ), where
    % R is the tangent sphere radius, and k is the aspheric factor:
    % 0 < k - oblate spheroid
    % k = 0 - sphere
    % -1 < k < 0 - prolate spheroid
    % k = -1 - parabola
    % k < -1 - hyperbola
    %
    % Member functions:
    %
    % l = Lens( r, D, R, k, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - diameter, 1x1 vector (outer) or 2x1 vector (inner, outer)
    % R - tangent sphere radius, [ Ry Rz ] vector for an astigmatic surface
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
        D = [ 0; 1 ]; % lens diameter (inner, outer)
     end
    
    methods
        function self = Lens( ar, aD, aR, ak, aglass )
            if nargin == 0
                return;
            end
            if size( aD, 1 ) < size( aD, 2 )
                aD = aD';
            end
            if size( aD, 1 ) == 1
                aD = [ 0; aD ];
            end
            self.r = ar;
            self.D = aD;
            self.R = aR;
            self.k = ak;
            self.glass = aglass;
            if ( self.D(2) / 2 / self.R(1) )^2 * ( 1 + self.k ) > 1
                warning('Lens Diameter is too large for its radius and apsheric parameter!' );
                self.D(1) = -1; % signal bad parameters
            end
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D(2) );
            if self.D(1) ~= 0
                fprintf( 'Inner diameter:\t %.3f\n', self.D(1) );
            end
            if length( self.R ) == 1
                fprintf( 'Curv. radius:\t %.3f\n', self.R );
            else
                fprintf( 'Y Curv. radius:\t %.3f\n', self.R(1) );
                fprintf( 'Z Curv. radius:\t %.3f\n', self.R(2) );
            end
            fprintf( 'Conic coefficient:\t %.3f\n', self.k );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self, color )
            % DISPLAY the lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            nrad = 50;
            rad = linspace( self.D(1) / 2, self.D(2) / 2, nrad );
            nang = 100;
            ang = linspace( 0, 2 * pi, nang );
            [ ang, rad ] = meshgrid( ang, rad );
            
            [ y, z ] = pol2cart( ang, rad );
            a = 1 + self.k;
            if length( self.R ) == 1
                r2yz = ( y.^2 + z.^2 ) / self.R^2;
                if a == 0 % paraboloid, special case
                    x = r2yz * self.R / 2;
                else
                    x = self.R * r2yz ./ ( 1 + real(sqrt( 1 - a * r2yz )) );
                end
            else % asymmetric conic
                r2yz = y.^2 / self.R(1) + z.^2 / self.R(2);
                if a == 0 % paraboloid, special case
                    x = r2yz / 2;
                else
                    x = r2yz ./ ( 1 + sqrt( 1 - a * ( y.^2 / self.R(1)^2 + self.R(2) / self.R(1) * z.^2 / self.R(2)^2 ) ) );
                end
            end
            

            S = [ x(:) y(:) z(:) ];
            
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
        
    end
    
end

