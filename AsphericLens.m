 classdef AsphericLens < GeneralLens
    % ASPHERICLENS Implements a lens surface given by a rotation of a conic curve
    % with additional polynomial terms, its sag given by
    % s = r^2 / R / ( 1 + sqrt( 1 - ( k + 1 ) * r^2 / R^2 ) ) + a(1) * r^2 +
    % a(2) * r^4 + ... + a(n) * r^2n, where
    % R is the tangent sphere radius, k is the aspheric factor:
    % 0 < k - oblate spheroid
    % k = 0 - sphere
    % -1 < k < 0 - prolate spheroid
    % k = -1 - parabola
    % k < -1 - hyperbola
    % and a(1) ... a(n) are the polynomial aspheric terms
    %
    % Member functions:
    %
    % l = AsphericLens( r, D, R, k, avec, glass ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - diameter
    % R - tangent sphere radius, [ Ry Rz ] vector for an astigmatic surface
    % k - conic coefficient, for astigmatic surface corresponds to the y-axis
    % avec - a vector of the polynomial terms
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
        avec = []; % vector of polynomial coeffeicient
    end
    
    methods
        function self = AsphericLens( ar, aD, aR, ak, aavec, aglass )
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
            if size( aavec, 1 ) > size( aavec, 2 )
                aavec = aavec';
            end
            self.avec = aavec;
            self.glass = aglass;
            self.funcs = 'asphlens';
            self.funch = str2func( self.funcs ); % construct function handle from the function name string
            if length( self.R ) == 1 % no astigmatism
                self.funca = [ self.R, self.R, self.k, self.avec ];
            else
                self.funca = [ self.R, self.k, self.avec ];
            end
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter:\t %.3f\n', self.D );
            fprintf( 'Curv. radius:\t %.3f\n', self.R );
            fprintf( 'Asphericity:\t %.3f\n', self.k );
            fprintf( 'Polynomial coefficients:' );
            disp( self.avec );
            fprintf( 'Material:\t %s | %s\n', self.glass{ 1 }, self.glass{ 2 } );
        end
        
        function h = draw( self, color )
            % DISPLAY the lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            h = self.draw@GeneralLens( color );
        end
        
        function rotate( self, rot_axis, rot_angle )
            self.rotate@Surface( rot_axis, rot_angle ); % rotate the surface members
            if abs( rot_angle ) > pi/2
                self.avec = -self.avec;
                self.funca = [ self.R, self.k, self.avec ];
            end
        end
    end
    
end

