classdef SurfaceGeneric < handle
    
    % SURFACEGENERIC defines universal arbitrary surfaces.
    % Surfaces can be described by a function (cartesian or radial) and
    % a surface topology array.
    % The surface normals are derived using a numerical gradient approach.
    
    % Member functions:
    %
    % b = a.copy - copy surface a to surface b
    %
    % rotate( rot_axis, rot_angle ) - rotate the surface by rot_angle
    % (radians) about the 1x3 rotation axis, inherited by derived classes.
    
    
    
    
    
    properties   % public access
        r = [ 0 0 0 ];  % location vector
        dim =[];        %dim vector (cart or polar) [outer inner]
        R = []; %legacy unused
        
        glass  = { 'vacuum' 'vacuum' }; % material in front and behind the surface
        
        %type: 'polar' or 'cartesian
        
        %shape description - cartesian f(x,y)
        shape_funch = [] % the corresponding function handle
        shape_funca = [] % argument list for the function
        
        %shape description - polar     f(theta, rho)
        shape_funch_polar = [] % the corresponding handle
        
        %additional surface profile term a[-1..1,-1..1]
        profile_array = []; %profile (function_handle or 2d array)
        
        %additional surface profile term a[-pi..pi,r]
        profile_array_polar = []; %profile (function_handle or 2d array)
        
        type = 'cart'; % 'polar','cart' -> impact on "draw/plot"
        
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
        
        function self = SurfaceGeneric( ar, aD, afunc, aglass, varargin)
            self.r = ar;
            self.dim = aD;
            self.shape_funch = afunc;
            self.glass = aglass;
            self.shape_funca = varargin;
            
        end
        
        
        
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
            
            if ~(rot_angle == 0)
                [self.rotax, self.rotang] =self.rodrigues_composition( self.rotax,self.rotang, rot_axis, rot_angle);
                self.n = rodrigues_rot( self.n, rot_axis, rot_angle );
            end
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
        
        function plot( self )
            y = linspace(-self.dim(1)/2,self.dim(1)/2,20);
            z = linspace(-self.dim(2)/2,self.dim(2)/2,20);
            [y,z]=meshgrid(y,z);
            [x,nrms] = self.eval(y,z);
            figure();
            subplot(2,2,1);
            surf(y,z,x);hold on;
            xlabel('y');ylabel('z');
            colorbar;
            title('total');
            
            subplot(2,2,2);
            surf(y,z,x);hold on;
            quiver3(y,z,x,nrms(:,:,2),nrms(:,:,3),nrms(:,:,1));
            axis equal; %normals look
            xlabel('y');ylabel('z');
            title('surface normals');
            
            subplot(2,2,3);
            [x] = self.eval_shape(y,z);
            surf(y,z,x);hold on;
            %quiver3(y,z,x,nrms1(:,:,2),nrms1(:,:,3),nrms1(:,:,1));
            title('shape f(y,z)');
            %axis equal; %view(90,-60);%normals look
            
            subplot(2,2,4);
            [x] = self.eval_profile(y,z);
            surf(y,z,x);hold on;
            %quiver3(y,z,x,nrms2(:,:,2),nrms2(:,:,3),nrms2(:,:,1));
            title('profile A(y,z)');
            %axis equal; %view(90,-60);%normals look v
            
            
        end
        
        
        
        function h = draw( self, color )
            % DISPLAY the lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            
            switch self.type
                case 'cart'
                    gridsize = 50;
                    y = linspace( -self.dim(1) / 2, self.dim(1) / 2, gridsize );
                    z = linspace( -self.dim(2) / 2, self.dim(2) / 2, gridsize );
                    [ y, z ] = meshgrid( y, z );
                case 'polar'
                    nrad = 50;
                    rad = linspace( 0, self.dim(1) / 2, nrad );
                    nang = 100;
                    ang = linspace( 0, 2 * pi, nang );
                    [ ang, rad ] = meshgrid( ang, rad );
                    [ y, z ] = pol2cart( ang, rad );
            end %switch
            
            x = self.eval( y, z );
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
        
        
        function profile_set(self, par)
            self.profile_array = par;
        end
        
        function [x]=eval_shape(self, yq,zq)
            x = -self.shape_funch(yq*2/self.dim(1),zq*2/self.dim(2));
            %remove imaginary numbers
             x(imag(x)~=0)=Inf;
        end
        %return profile for unit -1...1
        function [x]=eval_profile_unit(self, yq,zq)
               if ~isempty(self.profile_array)
                    ny = size(self.profile_array,2);
                    nz = size(self.profile_array,1);
                    yt=linspace(-1,1,ny);
                    zt=linspace(-1,1,nz);
                    x= interp2(yt,zt,self.profile_array,yq,zq,'linear',0);
               else
                    x=zeros(size(yq));
               end
            
        end
        
        function [x]=eval_profile(self, yq,zq)
            if ~isempty(self.profile_array)
                ny = size(self.profile_array,2);
                nz = size(self.profile_array,1);
                yt=linspace(-1,1,ny);
                zt=linspace(-1,1,nz);
                x= interp2(yt,zt,self.profile_array,yq*2/self.dim(1),zq*2/self.dim(2),'linear',0);         
            else
                x=zeros(size(yq));
            end
        end
        
        function [x,nrms] = eval(self,yq,zq)
            if nargout>1 %only calculated when asked
                [x1]=eval_shape(self, yq,zq);
                [x2]=eval_profile(self, yq,zq);
                x=x1+x2;       

                %grid for generating 2d surface for gradient (nrms)
                gridsize = 1e3;
                yt=linspace(-1,1,gridsize);
                zt=linspace(-1,1,gridsize);
                [yyt,zzt]= meshgrid(yt,zt);
                x_grid = self.shape_funch(yyt,zzt) + ...
                     self.eval_profile_unit(yyt,zzt) ;
                x_grid(imag(x_grid)~=0)=Inf; 
                 
                [gy,gz] = gradient(x_grid,2/(gridsize-1));
                %Remark: The gradient depends on correct scaling with
                %regards to y/z. When use a factor is needed.
                nrms_y = -interp2(yt,zt,gy,yq*2/self.dim(1),zq*2/self.dim(2),'cubic',0)./(self.dim(1)/2);
                nrms_z = -interp2(yt,zt,gz,yq*2/self.dim(1),zq*2/self.dim(2),'cubic',0)./(self.dim(2)/2);
                nrms = squeeze(cat(3,ones(size(nrms_y)),nrms_y,nrms_z));
                
               
                nrms = nrms./vecnorm(nrms,2,ndims(nrms));

            
            else
                [x1]=eval_shape(self, yq,zq);
                [x2]=eval_profile(self, yq,zq);
                x=x1+x2;
            end
            
        end
        
        
    end
    
end