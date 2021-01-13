classdef Surface < handle
    % SURFACE base class for optic elements
    %
    % Member functions:
    %
    % b = a.copy - copy surface a to surface b
    %
    % rotate( rot_axis, rot_angle ) - rotate the surface by rot_angle
    % (radians) about the 1x3 rotation axis, inherited by derived classes.
    %
    % Copyright: Yury Petrov, 2016
    %
    
    properties   % public access
        r = [ 0 0 0 ];  % location vector
        R = [];       % radius of the tangent sphere
        k = [];         % conic constant
        glass  = { 'vacuum' 'vacuum' }; % material in front and behind the surface
    end
    
    properties ( SetAccess = private )
        rotax = [ 1 0 0 ];  % rotation axis for rotation transformation
        rotang = 0;        % rotation angle about the rotax, radians      
        n = [ 1 0 0 ];      % orientation (normal) vector, can be set by a rotate function only
		
		

		
    end
    
	methods (Static)
		function [v_c, l_c] =rodrigues_composition( v_a,l_a, v_b, l_b )
			%combines two euler rodrigues rotations v_a*l_a and v_b*l_b 
			
		    l_c = 2*acos(  cos(l_a/2)*cos(l_b/2)-sin(l_a/2)*sin(l_b/2)*dot(v_a,v_b));
			v_c  = (sin(l_a/2)*cos(l_b/2)*v_a + cos(l_a/2)*sin(l_b/2)*v_b + ...
						 sin(l_a/2)*sin(l_b/2)*cross(v_a,v_b))./ sin(l_c/2);
		 end
	end
	
	
    methods
        
        function b = copy( self )
            b = feval( class( self ) );
            p = properties( self );
            for i = 1:length( p )
                b.( p{ i } ) = self.( p{ i } );
            end
        end

        function rotate( self, rot_axis, rot_angle )
            if abs( rot_angle ) > pi
                error( 'Rotation angle should be [ -pi pi ]!' );
            end
            % rotate the normal about the rot_axis by rot_angle (radians)
			
			

			%append rotation to existing rotation, Rodrigues rotation
			
			
			[self.rotax, self.rotang] =self.rodrigues_composition( self.rotax,self.rotang, rot_axis, rot_angle);
			
            self.n = rodrigues_rot( self.n, rot_axis, rot_angle );
        end
		
        
        
        function rotate_flip( self, rot_axis, rot_angle )
            if abs( rot_angle ) > pi
                error( 'Rotation angle should be [ -pi pi ]!' );
            end
            % rotate the normal about the rot_axis by rot_angle (radians)
            self.rotax = rot_axis;
            self.rotang = self.rotang + rot_angle;
            self.n = rodrigues_rot( self.n, rot_axis, rot_angle );
            if abs( self.rotang ) > pi/2 && ...
                    ~(isa(self,'Screen') || isa(self,'ScreenWavefront'))
                self.rotang = self.rotang - sign( self.rotang ) * pi;
                self.R = -self.R;  % surface upside-down
                self.glass = fliplr( self.glass ); % change material ordering
                self.n = -self.n; % make the surface normal point along the rays
            end
        end
 
    end
    
end