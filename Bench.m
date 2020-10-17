classdef Bench < handle
    % Bench class implements a system of optical elements
    % A complex optical system can be stored and manipulated as a whole by
    % making it a Bench instance.
    %
    % Member functions:
    %
    % b = Bench( obj )  - constructor function
    % INPUT:
    %   obj - an optical element, cell array of elements, or another bench
    % OUTPUT:
    %   b - bench object
    %
    % b.display() - displays bench b's information
    %
    % b.draw( rays, draw_fl, alpha, scale, new_figure_fl ) - draws bench b in the current axes
    % INPUT:
    %   rays - array of rays objects comprising full or partial light path
    %   draw_fl - display rays as 'arrows' (default), 'lines', or 'rays'
    %   alpha - opacity of optical surfaces from 0 to 1, default .33
    %   scale - scale of the arrow heads for the 'arrows' draw_fl
    %   new_figure_fl - 0, do not open, or 1, open (default)
    % 
    % a = b.copy() - copies bench b to bench a
    %
    % b.append( a, n ) - appends element a to bench b n times. n > 1
    % corresponds to multiple possible interactions (internal reflections
    % and such).
    %
    % b.prepend( a, n ) - prepends element a to bench b n times
    %
    % b.replace( ind, a ) - replaces an element with index ind on bench b with element a
    %
    % b.remove( inds ) - removes elements located at inds on bench b
    %
    % b.rotate( rot_axis, rot_angle, rot_fl ) - rotate the bench b with all its elements
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    %   rot_fl - (0, default) rotation of the bench elements wrt to the
    %   global origin, (1) rotation wrt to the bench geometric center
    %   
    % b.translate( tr_vec ) - translate the bench b with all its elements
    % INPUT:
    %   tr_vec - 1x3 translation vector
    %
    % rays_through = b.trace( rays_in, out_fl ) - trace rays through optical elements
    % on the bench b
    % INPUT:
    %   rays_in - incoming rays, e.g., created by the Rays() function
    %   out_fl  - 0 include even rays that missed some elements on the
    %   bench,  - 1 (default) exlude such rays
    % OUTPUT:
    %   rays_through - a cell array of refracted/reflected rays of the same
    %   length as the number of optical elements on the bench.
    %
    % Copyright: Yury Petrov, 2016
    %
  
    properties
        elem = {};     % cell array of optical elements
        cnt = 0;       % counter of elements in the system
    end
    
    methods
        function self = Bench( obj )
            % b = Bench( obj )  - constructor function
            % INPUT:
            %   obj - an optical element, cell array of elements, or another bench
            % OUTPUT:
            %   b - bench object
            if nargin == 0
                return;
            end
            
            if isa( obj, 'Bench' ) % if another Bench
                obj = obj.elem;    % extract elements
            end
            
            % append object(s) to the optical system
            nobj = length( obj );
            for i = 1 : nobj
                self.cnt = self.cnt + 1;
                if nobj == 1
                    self.elem{ self.cnt } = obj;
                elseif iscell( obj )   % other benches or cell arrays of Surfaces
                    self.elem{ self.cnt } = obj{ i };
                elseif isvector( obj ) % Rays
                    self.elem{ self.cnt } = obj( i );
                end
            end
        end
         
        function display( self )
            % b.display() - displays bench b's information
            for i = 1 : self.cnt
                obj = self.elem{ i };
                fprintf( '\n%s:\n', class( obj ) );
                obj.display;
            end
         end
        
        function draw( self, rays, draw_fl, alpha, scale, new_figure_fl )
            % b.draw( rays, draw_fl, alpha, scale, new_figure_fl ) - draws bench b in the current axes
            % INPUT:
            %   rays - array of rays objects comprising full or partial light path
            %   draw_fl - display rays as 'arrows' (default), 'lines', or 'rays'
            %   alpha - opacity of optical surfaces from 0 to 1, default .33
            %   scale - scale of the arrow heads for the 'arrows' draw_fl
            %   new_figure_fl - 0, do not open, or 1, open (default)
            if nargin < 6 || isempty( new_figure_fl )
                new_figure_fl = 1; % open a new figure by default
            end
            if nargin < 5 || isempty( scale )
                if nargin > 1
                    scale = ones( 1, length( rays ) );
                else
                    scale = 1;
                end
            else
                if length( scale ) == 1
                    scale = repmat( scale, 1, length( rays ) ); % make all ones
                elseif length( scale ) < length( rays )
                    if size( scale, 1 ) > size( scale, 2 )
                        scale = scale';
                    end
                    scale = [ scale ones( 1, length( rays ) - length( scale ) ) ]; % append ones 
                end
            end
            if nargin < 4 || isempty( alpha )
                alpha = 0.33;
            end
            if nargin < 3 || isempty( draw_fl )
                draw_fl = 'arrows';
            end
            if nargin < 2 || isempty( rays )
                rays = [];
            end
            
            if new_figure_fl == 1
                fname = dbstack;  % get debugging info
                
                if length(fname)>1
                    [ ~, fname ] = fname.name; % get the second (original) call function name
                else
                    fname ='';
                end
                
                
                figure( 'Name', [ 'OPTOMETRIKA: ' fname ], 'NumberTitle', 'Off', ...
                    'Position', [ 0 0 1024 1024 ], ...
                    'Color', 'k' );
            end
            hold on;
            for i = 1 : self.cnt
                obj = self.elem{ i };
                if isprop( obj, 'glass' ) && ( strcmp( obj.glass{1}, 'soot' ) || strcmp( obj.glass{2}, 'soot' ) )
                    color = [ .25 .25 .25 1 ];
                    obj.draw( color );
                else
                    color = [ 1 1 1 alpha ];
                    obj.draw( color );
                end
            end
            
            if ~isempty( rays )
                if strcmp( draw_fl, 'lines' ) || strcmp( draw_fl, 'clines' ) || strcmp( draw_fl, 'rays' ) % draw ray bundles as lines
                    if strcmp( draw_fl, 'lines' )
                        sym = '-';
                    else
                        sym = '*:';
                    end
                    for i = 1 : length( rays ) - 1
                        vis = ( rays( i ).I ~= 0 ) & ...
                                isfinite( sum( rays( i ).r.^2, 2 ) ) & ...
                                isfinite( sum( rays( i + 1 ).r.^2, 2 ) );  % visible rays
                        real = dot( rays( i + 1 ).r - rays( i ).r, rays( i ).n, 2 ) > 0; % real rays (vs. virtual for virtual image)
                        [ unique_colors, ~, ic ] = unique( rays( i ).color, 'rows' );
                        for j = 1 : size( unique_colors, 1 )
                            cvis = vis & real & ( ic == j );
                            plot3( [ rays( i ).r( cvis, 1 )';  rays( i + 1 ).r( cvis, 1 )' ], ...
                                   [ rays( i ).r( cvis, 2 )';  rays( i + 1 ).r( cvis, 2 )' ], ...
                                   [ rays( i ).r( cvis, 3 )';  rays( i + 1 ).r( cvis, 3 )' ], sym, 'Color', unique_colors( j, : ) );
                        end
                   end
                elseif strcmp( draw_fl, 'arrows' )
                    for i = 1 : length( rays )
                        rays( i ).draw( scale( i ) );
                    end
                end
            end
            
            %if new_figure_fl == 1
                gca.Clipping = 'off';
                axis equal vis3d off;
                %grid on;
                camlight( 'left' );
                camlight( 'right' );
                camlight( 'headlight' );
                view( -54, 54 );
                lighting phong;
                rotate3d on;
            %end
        end
        
        function b = copy( self )
            % a = b.copy() - copies bench b to bench a
            b = feval( class( self ) );
            b.cnt = self.cnt;
            for i = 1 : length( self.elem )
                b.elem{ i } = self.elem{ i }.copy;
            end
        end
        
        function append( self, obj, mult )
            % b.append( a, n ) - appends element a to bench b n times. n > 1
            % corresponds to multiple possible interactions (internal reflections
            % and such).
            if nargin < 3
                mult = 1;
            end
            if isa( obj, 'Bench' ) % if another Bench
                obj = obj.elem;    % extract elements
            end           
            % append object(s) to the optical system
            nobj = length( obj );
            for m = 1 : mult
                for i = 1 : nobj
                    self.cnt = self.cnt + 1;
                    if nobj == 1
                        self.elem{ self.cnt } = obj;
                    elseif iscell( obj )   % other benches or cell arrays of Surfaces
                        self.elem{ self.cnt } = obj{ i };
                    elseif isvector( obj ) % Rays
                        self.elem{ self.cnt } = obj( i );
                    end
                end
            end
        end
        
        function prepend( self, obj, mult )
            % b.prepend( a, n ) - prepends element a to bench b n times
            if nargin < 3
                mult = 1;
            end
            if isa( obj, 'Bench' ) % if another Bench
                obj = obj.elem;    % extract elements
            end         
            self.elem = fliplr( self.elem ); % reverse element direction temporarily
            % prepend object(s) to the optical system
            nobj = length( obj );
            for m = 1 : mult
                for i = nobj : -1 : 1 % append in the opposite order
                    self.cnt = self.cnt + 1;
                    if nobj == 1
                        self.elem{ self.cnt } = obj;
                    elseif iscell( obj )   % other benches or cell arrays of Surfaces
                        self.elem{ self.cnt } = obj{ i };
                    elseif isvector( obj ) % Rays
                        self.elem{ self.cnt } = obj( i );
                    end
                end
            end
            self.elem = fliplr( self.elem ); % restitute the original order
        end
        
        function replace( self, ind, obj )
            % b.replace( ind, a ) - replaces an element with index ind on bench b with element a
            self.elem{ ind } = obj;
        end
        
         function remove( self, inds )
             % b.remove( inds ) - removes elements located at inds on bench b
             if self.cnt == 0
                 error( 'The bench is already empty!' );
             else
                self.elem( inds ) = [];
                self.cnt = self.cnt - length( inds );
             end
         end
         
         function rotate( self, rot_axis, rot_angle, rot_fl )
             % b.rotate( rot_axis, rot_angle, rot_fl ) - rotate the bench b with all its elements
             % INPUT:
             %   rot_axis - 1x3 vector defining the rotation axis
             %   rot_angle - rotation angle (radians)
             %   rot_fl - (0, default) rotation of the bench elements wrt to the
             %   global origin, (1) rotation wrt to the bench geometric center
             if nargin < 4
                 rot_fl = 0;
             end
             cntr = [ 0 0 0 ];
             if rot_fl == 1 % rotate around the geometric center of the bench
                 for i = 1 : self.cnt % loop through the optic system
                     cntr = cntr + self.elem{ i }.r;
                 end
                 cntr = cntr / self.cnt;
             end
            % rotate bench elements
            for i = 1 : self.cnt % loop through the optic system
                self.elem{ i }.rotate( rot_axis, rot_angle ); % rotate normal
                self.elem{ i }.r = cntr + rodrigues_rot( self.elem{ i }.r - cntr, rot_axis, rot_angle ); % rotate position
            end
            if abs( rot_angle ) > pi/2 % reverse order in which the elements are encountered by rays
                self.elem = fliplr( self.elem );
            end
        end

        function translate( self, tr_vec )
            % b.translate( tr_vec ) - translate the bench b with all its elements
            % INPUT:
            %   tr_vec - 1x3 translation vector
            for i = 1 : self.cnt % loop through the optic system
                self.elem{ i }.r = self.elem{ i }.r + tr_vec; % translate position
            end
        end
       
        function rays = trace( self, rays_in, out_fl )
            % rays_through = b.trace( rays_in, out_fl ) - trace rays through optical elements
            % on the bench b
            % INPUT:
            %   rays_in - incoming rays, e.g., created by the Rays() function
            %   out_fl  - 0 include even rays that missed some elements on the
            %   bench,  - 1 (default) exlude such rays
            % OUTPUT:
            %   rays_through - a cell array of refracted/reflected rays of the same
            %   length as the number of optical elements on the bench.
            if nargin < 3
                out_fl = 1; % exclude rays which miss elements of the bench
            end
            rays( 1, self.cnt + 1 ) = Rays; % allocate cnt instances of Rays
            rays( 1 ) = rays_in;
            for i = 1 : self.cnt % loop through the optic system
                rays( i + 1 ) = rays( i ).interaction( self.elem{ i }, out_fl );
                
%                 if(0) %helpful debug output of traced rays
%                     disp(['Tracing through element ' num2str(i) ': ' class(self.elem{i}) ]);
%                     %.r -> next origin
%                     %.n -> next direction
%                     disp(['r: ' num2str(mean(rays(i).r,'omitnan')) '->' num2str(mean(rays(i+1).r,'omitnan'))]) 
%                     disp(['n: ' num2str(mean(rays(i).n,'omitnan')) '->' num2str(mean(rays(i+1).n,'omitnan'))]) 
%                 end         
                
            end
        end % trace()
 
function rays = trace_recursive( self, rays_in, out_fl )
            % rays_through = b.trace( rays_in, out_fl ) - trace rays through optical elements
            % on the bench recursively, finding hitting points
            % INPUT:
            %   rays_in - incoming rays, e.g., created by the Rays() function
            %   out_fl  - 0 include even rays that missed some elements on the
            %   bench,  - 1 (default) exlude such rays
            % OUTPUT:
            %   rays_through - a cell array of refracted/reflected rays of the same
            %   length as the number of optical elements on the bench.
            if nargin < 3
                out_fl = 1; % exclude rays which miss elements of the bench
            end
            
            warning('off', 'OptoMetrika:refrIndexMismatch');
            
            %We do not know yet how many steps it will take
            %rays( 1, self.cnt + 1 ) = Rays; % allocate cnt instances of Rays
            
            
            rays( 1 ) = rays_in;
            
            jj = 1;
            
            abort = false;

            while(~abort)
            
                %Calculate possible intersection with all elements
                for i = 1 : self.cnt % loop through the optic system
                    tmp_rays( i ) = rays( jj ).interaction( self.elem{ i }, out_fl );  
                end

                %tmp_rays( i ) = rays( jj-1 ).interaction( self.elem{ 3 }, out_fl );

                %Choose the shortest & valid one
                gpl=[tmp_rays.gpl];
                
                %is not valid ray that has not missed
                gpl(~([tmp_rays.I] >=eps)) = NaN;
                
                %avoid the current element(distance = 0)
                % and elements in negative direction (distance < 0)
                gpl(gpl <0) = NaN;
                
                [~,indx]=min(gpl,[],2);

                %Abort condition. Check how many valid ray's were in the last
                valid_rays_cnt=sum(([gpl] >=0),'all');
                abort = (valid_rays_cnt==0);
                
                if(~abort)
                              
                    %Build the new ray vector depending on selected index
                    tmp_rays_new=Rays();
                    for i = 1 : self.cnt
                        sel = find(indx==i);                  
                        tmp_rays_new=tmp_rays_new.append(tmp_rays(i).subset(sel));
                    end 

                    jj=jj+1;
                    rays(jj) = tmp_rays_new;
                
                end %~abort
                
%                 if(true) %helpful debug output of traced rays
%                     %.r -> next origin
%                     %.n -> next direction
%                     disp(['r: ' num2str(mean(rays(jj-1).r,'omitnan')) '->' num2str(mean(rays(jj).r,'omitnan'))]) 
%                     disp(['n: ' num2str(mean(rays(jj-1).n,'omitnan')) '->' num2str(mean(rays(jj).n,'omitnan'))]) 
%                 end      
            end% while ~abort
     
            warning('on', 'OptoMetrika:refrIndexMismatch');
        end % trace_recursive()
        
        
    end
    
end

