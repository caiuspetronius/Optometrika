classdef GeneralLens < Surface
    % GENERALLENS Implements a lens surface of an arbitrary shape. This
    % class requires iterative search of intersections with a ray,
    % therefore it works slower than Lens. 
    %
    % Member functions:
    %
    % l = GeneralLens( r, D, func, glass, varargin ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - diameter
    % func - function name string
    % glass - 1 x 2 cell array of strings, e.g., { 'air' 'acrylic' }
    % varargin - an arbitrary number of parameters required by func.
    % OUTPUT:
    % l - lens surface object
    %
    % display() - displays the object's information
    %
    % draw() - draws the object in the current axes
    %
    % rotate( rot_axis, rot_angle ) - rotate the surface by rot_angle
    % (radians) about the 1x3 rotation axis.
    % 
    % Copyright: Yury Petrov, 2016
    %
    
    properties
        D = [ 0; 1 ]      % lens diameter (inner, outer)
        funcs = '' % lens surface function name string
        funch = [] % the corresponding function handle
        funca = [] % argument list for the function
     end
    
    methods
        function self = GeneralLens( ar, aD, afunc, aglass, varargin )
            if nargin == 0
                return;
            end
            if size( aD, 1 ) < size( aD, 2 )
                aD = aD';
            end
            if size( aD, 1 ) == 1
                aD = [ 0; aD ]; % assume inner radius = 0
            end
            self.r = ar;
            self.D = aD;
            self.funcs = afunc;
            self.funch = str2func( afunc ); % construct function handle from the function name string
            self.glass = aglass;
            self.funca = varargin;
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D(2) );
            if self.D(1) ~= 0
                fprintf( 'Inner diameter:\t %.3f\n', self.D(1) );
            end
            fprintf( 'Surface function:\t %s\n', self.func );
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
            x = self.funch( y, z, self.funca, 0 );
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

