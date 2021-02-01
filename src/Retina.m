classdef Retina < Surface
    % RETINA implements a spherical/ellipsoidal screen
    %
    % Member functions:
    %
    % p = Retina( r, R, k ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % R - tangent sphere radius
    % k - conic coefficient for the surface
    % OUTPUT:
    % p - retina object
    %
    % p.display() - displays the retina p information
    %
    % p.draw() - draws the retina p in the current axes
    % 
    % p.rotate( rot_axis, rot_angle ) - rotate the retina
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    %
    % Copyright: Yury Petrov, 2016
    %
  
    properties
        D = 1;   % retina diameter
        ang = pi/2;  % aperture angle
        azbins = 512; % number of azimuth bins
        elbins = 512; % number of elevation bins
        image = [];
    end
    
    properties ( SetAccess = private )
    end
    
    methods
        function self = Retina( ar, aR, ak, aang, aazbins, aelbins )
            if nargin == 0
                return;
            end
            if nargin < 6
                aelbins = 512;
            end
            if nargin < 5
                aazbins = 512;
            end
            if nargin < 4
                aang = 0.69 * pi/2; % to fit the lens
            end
            if nargin < 3
                ak = 0;
            end
            if nargin < 2
                error( 'At least position and radius Retina  parameters must be specified!' );
            end
            self.r = ar;
            self.R = aR;
            self.k = ak;
            self.ang = aang;
            self.azbins = aazbins;
            self.elbins = aelbins;
            if ak <= -1 % not a sphere or ellipsoid
               error( 'Aspheric parameter has to be larger than -1 for Retina' );
            end
            self.D = 2 * abs( self.R ) ./ sqrt( 1 + self.k );
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D );
            fprintf( 'Curv. radius:\t %.3f\n', self.R );
            fprintf( 'Asphericity:\t %.3f\n', self.k );
        end
        
        function h = draw( self, color )
            % DISPLAY the spherical surface
            if nargin < 2
                color = [ 0.2 0.2 0.2 1 ]; % dull gray opaque color
            end
            theta = linspace( -pi, pi, self.azbins );     % azimuth
            phi   = linspace( -pi/2, self.ang, self.elbins )'; % elevation
            cosphi = cos( phi ); cosphi(1) = 0;
            sintheta = sin( theta ); sintheta(1) = 0; sintheta( end ) = 0;
            y = cosphi * cos( theta );
            z = cosphi * sintheta;
            x = sin( phi ) * ones( 1, length( theta ) );
            S = self.R * [ ( x(:) + 1 ) / ( 1 + self.k ), y(:) / sqrt( 1 + self.k ), z(:) / sqrt( 1 + self.k ) ]; % apex at the orgin
            
            % rotate and shift
            if self.rotang ~= 0
                S = rodrigues_rot( S, self.rotax, self.rotang );
            end
            x(:) = S( :, 1 ) + self.r( 1 );
            y(:) = S( :, 2 ) + self.r( 2 );
            z(:) = S( :, 3 ) + self.r( 3 );
            
            if ~isempty( self.image )
                c = self.image;
            else
                c = repmat( reshape( color( 1:3 ), [ 1 1 3 ] ), size( x, 1 ), size( x, 2 ), 1 );
            end
            h = surf( x, y, z, c, ...
                'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', 'FaceAlpha', color(4), ...
                'AmbientStrength', 1, 'SpecularStrength', 0 );
            colormap summer;
        end

    end
    
end

