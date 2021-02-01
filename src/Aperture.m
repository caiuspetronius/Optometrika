classdef Aperture < Surface
    % APERTURE defines a circular or rectangular opening
    %
    % Member functions:
    %
    % a = Aperture( r, D ) - object constructor
    % INPUT:
    % r - 1x3 position vector
    % D - 2x1 vector (inner diameter, outer diameter) or 4x1 vector (inner
    % w, innher h, outer w, outer h)
    % OUTPUT:
    % a - aperture object
    %
    % a.display() - displays the aperture a's information
    %
    % draw() - draws the aperture a in the current axes
    % 
    % a.rotate( rot_axis, rot_angle ) - rotate the aperture
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    % 
    % Copyright: Yury Petrov, 2016
    %
    
    properties
        D = [ 1; 2 ]
    end
    
    methods
        function self = Aperture( ar, aD )
            self.glass = { 'air' 'soot' };
            if nargin == 0
                return;
            end
            self.r = ar;
            if size( aD, 1 ) < size( aD, 2 )
                aD = aD';
            end
            self.D = aD;
            if size( self.D, 1 ) ~= 2 && size( self.D, 1 ) ~= 4 
                error( 'Aperture dimensions have to be a 2x1 or 4x1 vector' );
            elseif size( self.D, 1 ) == 2 && self.D(2) < self.D(1)
                error( 'Outer radius has to be larger than the inner radius' );
            elseif size( self.D, 1 ) == 4 && ( self.D(1) > self.D(3) || self.D(2) > self.D(4) )
                error( 'Outer aperture dimensions have to be larger than inner dimensions' );
            end
        end
        
        function display( self )
            fprintf( 'Position:\t [%.3f %.3f %.3f]\n', self.r );
            fprintf( 'Orientation:\t [%.3f %.3f %.3f]\n', self.n );
            fprintf( 'Diameter in:\t %.3f\n', self.D(1) );
            fprintf( 'Diameter out:\t %.3f\n', self.D(2) );
        end
        
        function h = draw( self, color )
            if size( self.D, 1 ) == 2 % circular aperture
                nrad = 2;
                rad = linspace( self.D(1) / 2, self.D(2) / 2, nrad );
                nang = 100;
                ang = linspace( 0, 2 * pi, nang );
                [ ang, rad ] = meshgrid( ang, rad );
                
                y = rad .* cos( ang );
                z = rad .* sin( ang );
                x = zeros( size( y ) );
                S = [ x(:) y(:) z(:) ];
                
                % rotate and shift
                if self.rotang ~= 0
                    S = rodrigues_rot( S, self.rotax, self.rotang );
                end
                x(:) = S( :, 1 ) + self.r( 1 );
                y(:) = S( :, 2 ) + self.r( 2 );
                z(:) = S( :, 3 ) + self.r( 3 );
                % draw
                c = 0.25 * ones( size( x, 1 ), size( x, 2 ), 3 );
                h = surf( x, y, z, c, 'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', 'interp', ...
                'AmbientStrength', 0.25, 'SpecularStrength', 0.25 );
            else % rectangular aperture
                w = self.D(1);
                h = self.D(2);
                W = self.D(3);
                H = self.D(4);
                V = [ 0 0 0
                      0 W 0
                      0 W H
                      0 0 H ];
                v = [ 0 0 0
                      0 w 0
                      0 w h
                      0 0 h ];
                v( :, 2 ) = v( :, 2 ) + ( W - w )/2;
                v( :, 3 ) = v( :, 3 ) + ( H - h )/2;
                v = [ V; v ];
                v( :, 1 ) = v( :, 1 ) + self.r( 1 );
                v( :, 2 ) = v( :, 2 ) + self.r( 2 ) - W/2;
                v( :, 3 ) = v( :, 3 ) + self.r( 3 ) - H/2;
                
                f = [ 1 2 6 5
                      2 3 7 6
                      3 4 8 7
                      4 1 5 8 ];
                  
                % rotate and shift
                S = [ v( :, 1 ) v( :, 2 ) v( :, 3 ) ];
                if self.rotang ~= 0
                    S = rodrigues_rot( S, self.rotax, self.rotang );
                end
                v = S + repmat( self.r, size( v, 1 ), 1 );
                  
                % draw
                h = patch( 'Faces', f, 'Vertices', v, 'EdgeColor', 'none', 'FaceLighting','phong', 'FaceColor', [ .25 .25 .25 ], ...
                      'AmbientStrength', 0.25, 'SpecularStrength', 0.25);
            end
        end
    end
    
end

