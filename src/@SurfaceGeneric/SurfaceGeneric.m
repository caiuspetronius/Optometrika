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
    
    
    % For 'spherical' the center of the sphere is at [0,0,0]
    % Definition: phi-azimuth- [-pi..pi],theta-elevation- [-pi/2..pi/2]
    
    
    
    properties   % public access
        r = [ 0 0 0 ];  % location vector
        dim =[];        %dim vector (cart or polar) [outer inner]
        R = []; 
        
        glass  = { 'vacuum' 'vacuum' }; % material in front and behind the surface
        
        %type: 'polar' or 'cartesian
        
        %shape description - cartesian f(x,y)
        shape_funch = [] % the corresponding function handle
        shape_funca = [] % argument list for the function
        
        %shape description - polar     f(theta, rho)
        shape_funch_polar = [] % the corresponding handle
        
        %additional surface profile term a[-1..1,-1..1]
         %A[y,z] for cartesian
        profile_array = []; %profile (function_handle or 2d array)
        
                             
        %additional surface profile term a[-pi..pi,r]
        profile_array_polar = []; %profile (function_handle or 2d array)
        %A[az,el] for spherical
             
        
        type = 'cart'; % 'polar','cart','spherical' defines what type of function
        %cart: f(x,y)
        %polar: f(theta, rho)
        %spherical: f(theta, phi, r)
        
        diffdist = 1e-5; %numeric differentiation distance for gradient
  
        
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
            
             %cartesian
             %upper and lower limits: -1 / 1
             y = linspace(-self.dim(1)/2,self.dim(1)/2,20);
             z = linspace(-self.dim(2)/2,self.dim(2)/2,20); 
             
            [y,z]=meshgrid(y,z);


            
            [x,nrms] = self.eval(y(:),z(:));
            x= reshape(x,size(y)); %back to mesh
            nrms= reshape(nrms,[size(y),3]); %back to mesh
            
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
        
       function plot_spherical( self )
             %spherical
             phi = linspace(-pi,pi,20);
             theta = linspace(-pi/2,pi/2,20); 
             
            [phi ,theta]=meshgrid(phi ,theta);
            [rr] = self.eval_spherical(phi ,theta);
            figure();
            subplot(2,2,1);
            surf(phi ,theta,rr);hold on;
            xlabel('phi(az),rad');ylabel('theta(el),rad');
            colorbar;
            title('total');

            
            subplot(2,2,3);
            [rr] = self.eval_shape_spherical(phi ,theta);
            surf(phi ,theta,rr);hold on;
            title('shape f(y,z)');
            xlabel('phi(az),rad');ylabel('theta(el),rad');
            subplot(2,2,4);
            [rr] = self.eval_profile_spherical(phi ,theta);
            surf(phi ,theta,rr);hold on;
            title('profile A(y,z)'); 
            xlabel('phi(az),rad');ylabel('theta(el),rad');
        end
        
        
        
        function h = draw( self, color )
            % DISPLAY the lens surface
            if nargin < 2
                color = [ 1 1 1 .5 ];
            end
            
            if ~isempty(self.R)
               self.type='spherical'; 
            end
            
            switch self.type
                case 'cart'
                    gridsize = 50;
                    y = linspace( -self.dim(1) / 2, self.dim(1) / 2, gridsize );
                    z = linspace( -self.dim(2) / 2, self.dim(2) / 2, gridsize );
                    [ y, z ] = meshgrid( y, z );

                    x = self.eval( y(:), z(:) );
                    x = reshape(x,size(y));
                case 'polar'
                    nrad = 50;
                    rad = linspace( 0, self.dim(1) / 2, nrad );
                    nang = 100;
                    ang = linspace( 0, 2 * pi, nang );
                    [ ang, rad ] = meshgrid( ang, rad );
                    [ y, z ] = pol2cart( ang, rad );
                    x = self.eval( y, z );
                case 'spherical'
                    nphi = 51;
                    ntheta = 51;
                    phi = linspace( -pi,pi, nphi );
                    theta = linspace( -pi/2,pi/2, ntheta );
                    [ phi, theta ] = meshgrid( phi, theta );
                    %rotate so that (0,0) is towards -x, the default
                    
                    %phi=0;theta=0;
                    rr = self.eval_spherical( phi, theta )+self.R;
                  
                    %[y,z,x] = sph2cart(phi,theta,rr);  
                    [x,y,z]=sph2cart(phi,theta,rr);  
                    %Rotate around z by 180? (face towards -x)
                    x=-x;
                    y=-y;           
                    
            end %switch
            
            
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
     
        
        function [x]=eval_shape_spherical(self, varphi_q,theta_q)
            x = self.shape_funch(theta_q,varphi_q);
            %remove imaginary numbers
            x(imag(x)~=0)=1e20;
        end
        
        function [x]=eval_profile_spherical(self,varphi_q,theta_q)
               if ~isempty(self.profile_array)
                    n_phi = size(self.profile_array,1); %az
                    n_theta = size(self.profile_array,2); %el
                    phi_t=linspace(-pi,pi,n_phi); %az
                    theta_t=linspace(-pi/2,pi/2,n_theta); %el
                    x= interp2(phi_t',theta_t',self.profile_array',varphi_q,theta_q,'linear',0);
               else
                    x=zeros(size(varphi_q));
               end
        end
        
        
        function [x]=eval_shape(self, yq,zq)
            x = -self.shape_funch(yq*2/self.dim(1),zq*2/self.dim(2));
            %remove imaginary numbers
             x(imag(x)~=0)=NaN;
        end
        %return profile for unit -1...1
        function [x]=eval_profile_unit(self, yq,zq)
               if ~isempty(self.profile_array)
                    ny = size(self.profile_array,2);
                    nz = size(self.profile_array,1);
                    yt=linspace(-1,1,ny);
                    zt=linspace(-1,1,nz);
                    x= interp2(yt,zt,self.profile_array,yq,zq,'makima');
                    %alternative: switch to griddedInterpolant for more
                    %control on intra/extrapolation
                    
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
                x= interp2(yt,zt,self.profile_array,yq*2/self.dim(1),zq*2/self.dim(2),'makima');         
            else
                x=zeros(size(yq));
            end
        end
        
        function [r] = eval_spherical(self, varphi_q,theta_q)    
            [r1]=eval_profile_spherical(self,varphi_q,theta_q);
            [r2]=eval_shape_spherical(self,varphi_q,theta_q);
            r=r1+r2;
        end
        
        
        function [x,nrms] = eval(self,yq,zq,sign_dirx)
            
            % cart to spherical
            if ~isempty(self.R)
                
                
                if ~exist('sign_dirx','var')
                    sign_dirx=1;
                end
                %when calculating the the spherical coordinate x
                %based on r there might be a slight error
                %for highly deviated surfaces
                %because it is performed for ideal sphere
                       
                %there are two solutions for intersection
                x1=sign_dirx.*sqrt(self.R.^2-yq.^2-zq.^2);
                x1(imag(x1)~=0)=NaN;
                
        
                
                theta= asin(zq./self.R); %el
                phi = atan2(yq,x1);      %az
                
                %theta= asin(x./self.R); %el
                %phi = atan2(yq,x);      %az
                
                phi(imag(phi)~=0)=NaN;
                theta(imag(theta)~=0)=NaN;
                
                [rr] = eval_spherical(self, phi,theta);
                %[x,~,~] = sph2cart(phi,theta,rr+self.R);
                
                [x,~,~] = sph2cart(phi,theta,rr+self.R);                   
                %Rotate around z by 180? (face towards -x)
                x=-x;
   
                
                
        
            else
                %cartesian
                [x1]=eval_shape(self, yq,zq);
                [x2]=eval_profile(self, yq,zq);
                x=x1+x2;
                
            end
            
               if nargout>1 %only calculated when asked
                %cartesian           
                 nrms =eval_nrms(self,yq,zq);           
               end
            
        end
        
        function [nrms] = eval_nrms(self,yq,zq)
                % calculate the normal vectors for 2 point pairs along x and y
                % normal vector is tangential using numerical
                % differentiation
                
                y1=yq+self.diffdist;
                y2=yq-self.diffdist;
                z1=zq+self.diffdist;
                z2=zq-self.diffdist;
                
                x= self.eval(yq,zq);
                x1=self.eval(y1,zq);
                x2=self.eval(y2,zq);
                x3=self.eval(yq,z1);
                x4=self.eval(yq,z2);
                
                v1=[x1,y1,zq];
                v2=[x2,y2,zq];
                v3=[x3,yq,z1];
                v4=[x4,yq,z2];

                vy = v2-v1;
                vz = v4-v3;
                nrms = cross(vz,vy);
                nrms = -nrms./vecnorm(nrms,2,2);
        end
        
        
    end
    
end