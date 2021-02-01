classdef Rays
    % RAYS Implements a ray bundle
    % Note that for easy copying Rays doesn't inherit from handle
    %
    % Member functions:
    %
    % r = Rays( n, geometry, r, dir, D, pattern, glass, wavelength, color, diopter ) - object constructor
    % INPUT:
    %   n - number of rays in the bundle
    %   geometry - For geometry 'collimated', r defines rays origins while dir - 
    % their direction. For geometry 'source', r defines position 
    % of the point source, and dir - direction along which rays propagate. For geometry
    % 'vergent' r defines rays orignis, dir - their average direction, while diopter 
    % defines the convergence/divergence of the rays in diopters.
    %   r - 1x3 bundle source position vector
    %   dir - 1x3 bundle direction vector
    %   D - diameter of the ray bundle (at distance 1 if geometry = 'source' )
    %   pattern - (optional) pattern of rays within the bundle: 'linear' = 'linearY', 'linearZ', 'hexagonal'
    %   'square', 'sphere' or 'random', hexagonal by default
    %   glass - (optional) material through which rays propagate, 'vacuum' by
    % default
    %   wavelength - (optional) wavelength of the ray bundle, meters, 557.7
    % nm by default
    %   color - (optional) 1 x 3 vector defining the color with which to draw
    % the ray, [ 0 1 0 ] (green) by default
    %
    % OUTPUT:
    %   r - ray bundle object
    %
    % r.draw( scale ) - draws the ray bundle r in the current axes as arrows
    % INPUT:
    %   scale - (optional) the arrow length, 1 by default
    %
    % [ rays_out, nrms ] = r.intersection( surf ) - finds intersection
    % of the ray bundle r with a surface
    % INPUT:
    %   surf - the surface
    % OUTPUT:
    %   rays_out - rays_out.r has the intersection points, the remaining
    % structure parameters might be wrong
    %   nrms - 1x3 normal vectors to the surface at the intersection points
    %
    % rays_out = r.interaction( surf ) - finishes forming the outcoming
    % ray bundle, calculates correct directions and intensities
    % INPUT:
    %   surf - surface
    % OUTPUT:
    %   rays_out - outcoming ray bundle
    %
    % r = r.append( r1 ) - appends to bundle r bundle r1
    % INPUT:
    %   r1 - the appended ray bundle
    % OUTPUT:
    %   r - the resulting ray bundle
    % 
    % sr = r.subset( inices ) - subset of rays in bundle r
    % INPUT:
    %   indices - subset's indices in the bundle r
    % OUTPUT:
    %   sr - the resulting new ray bundle
    % 
    % r = r.truncate() - truncate all rays with zero intensity from the bundle r.
    %
    % [ av, dv, nrays ] = r.stat() - return statistics on the bundle r
    % OUTPUT:
    %   av - 1x3 vector of the mean bundle position
    %   dv - standard deviation of the ray positions in the bundle
    %   nrays - number of rays with non-zero intensity in the bundle
    %
    % [ x0, cv, ax, ang, nrays ] = r.stat_ellipse() - fit a circumscribing ellipse
    % to the bundle r in the YZ plane
    % OUTPUT:
    %   x0 - 1x3 vector of the ellipse center
    %   cv - bundle covariance matrix
    %   ax - 1x2 vector of the ellipse half-axes lengths
    %   ang - angle of rotation of the ellipse from the longer axis being
    %   oriented along the Y axis.
    %   nrays - number of rays with non-zero intensity in the bundle
    %
    % r2 = dist2rays( p ) - returns squared distances from point p to all rays
    % INPUT:
    %   p - 1x3 vector
    % OUTPUT:
    %   r2 - nrays x 1 vector of squared distances
    %
    % [ f, ff ] = r.focal_point() - find a focal point of the bundle. The focal
    % point f is defined as the mean convergence distance of the bundle
    % rays. ff gives the residual bundle crossection (intensity weighted std).
    % OUTPUT:
    %    f - 1x3 vector for the focal point 
    %
    % Copyright: Yury Petrov, 2016
    %
    
    
    % Flag Function to change to higher precision numeric solution
    % of generic lenses. 
    methods (Static)   
      function ret = setget_solver_flag(par)
         persistent solver_flag;        
         if nargin
            solver_flag = par;
         elseif(isempty(solver_flag))
            solver_flag='default'; %default
         end
         ret = solver_flag;
      end
   end
    
   
    
    properties
        r = []; % a matrix of ray starting positions
        n = []; % a matrix of ray directions
        w = [];  % a vector of ray wavelengths
        I = [];         % a vector of ray intensities
        nrefr = [];    % a vector of current refractive indices
        att = [];       % a vector of ray attenuations
        color = [];   % color to draw the bundle rays
        cnt = 0;      % number of rays in the bundle
        opl = [];  % a vector of summed traced ray optical path length
        gpl = [];  % a vector of last traced ray geometrical path length
    end
    
    methods            
        function self = Rays( cnt, geometry, pos, dir, diameter, rflag, glass, wavelength, acolor, adiopter ) % constructor of ray bundles
            % Constructs a ray bundle comprising 'cnt' rays. For geometry 
            % 'collimated', 'pos' defines rays origins, while 'dir' - 
            % their direction. For geometry 'source', 'pos' defines position 
            % of the point source, 'dir' - direction along which rays form a 
            % linear, hexagonal, square, or random pattern (specified by 'rflag') of 
            % the size specified by 'diameter' at distance 1. For geometry
            % 'vergent' 'pos' defines rays orignis, 'dir' - their average
            % direction, while 'adiopter' defines the convergence/divergence
            % of the rays in diopters.
            
            if nargin == 0 % used to allocate arrays of Rays
                return;
            end
            if nargin < 10 || isempty( adiopter )
               diopter = 0;
            else
               diopter = adiopter;
            end
            if nargin < 9 || isempty( acolor )
                self.color = [ 0 1 0 ];
            else
                self.color = acolor;
            end
            if nargin < 8 || isempty( wavelength )
                self.w = 5300e-10; % green by default
            else
                self.w = wavelength;
            end
            if nargin < 7 || isempty( glass )
                glass = 'vacuum';
            end
            if nargin < 6 || isempty( rflag )
                rflag = 'hexagonal'; % hexagonal lattice of rays
            end
            if nargin < 5 || isempty( diameter )
                diameter = 1;
            end
            if nargin < 4 || isempty( dir )
                dir = [ 1 0 0 ];
            end
            if nargin < 3 || isempty( pos )
                pos = [ 0 0 0 ];
            end
            if nargin < 2 || isempty( geometry )
                geometry = 'collimated';
            end
                    
            % normalize direction and rotate positions to the plane
            % orthogonal to their direction
            dir = dir ./ norm( dir );
            ex = [ 1 0 0 ];
            ax = cross( ex, dir );

            if strcmp( rflag, 'linear' ) || strcmp( rflag, 'linearY' ) % extend along y-axis
                p( :, 1 ) = linspace( -diameter/2, diameter/2, cnt ); % all rays starting from the center
                p( :, 2 ) = 0;
            elseif strcmp( rflag, 'linearZ' ) % extend along z-axis
                p( :, 1 ) = zeros( 1, cnt );
                p( :, 2 ) = linspace( -diameter/2, diameter/2, cnt ); % all rays starting from the center
            elseif strcmp( rflag, 'random' )
                cnt1 = round( cnt * 4 / pi );
                p( :, 1 ) = diameter * ( rand( cnt1, 1 ) - 0.5 ); % horizontal positions
                p( :, 2 ) = diameter * ( rand( cnt1, 1 ) - 0.5 ); % vertical positions
                p( p( :, 1 ).^2 + p( :, 2 ).^2 > diameter^2 / 4, : ) = []; % leave rays only within the diameter
            elseif strcmp( rflag, 'hexagonal' )
                % find the closest hexagonal number to cnt
                cnt1 = round( cnt * 2 * sqrt(3) / pi );
                tmp = (-3 + sqrt( 9 - 12 * ( 1 - cnt1 ) ) ) / 6;
                cn( 1 ) = floor( tmp );
                cn( 2 ) = ceil(  tmp );
                totn = 1 + 3 * cn .* ( 1 + cn );
                [ ~, i ] = min( abs( totn - cnt1 ) );
                cn = cn( i );
                % generate hexagonal grid
                p = [];
                for i = cn : -1 : -cn % loop over rows starting from the top
                    nr = 2 * cn + 1 - abs( i ); % number in a row
                    hn = floor( nr / 2 );
                    if rem( nr, 2 ) == 1
                        x = ( -hn : hn )';
                    else
                        x = ( -hn : hn - 1 )' + 1/2;
                    end
                    p = [ p; [ x, i * sqrt( 3 ) / 2 * ones( nr, 1 ) ] ]; % add new pin locations
                end
                if cn > 0
                    p = p * diameter / 2 / cn * 2 / sqrt( 3 ); % circubscribe the hexagon by an inward circle
                end
%                 if cn > 2 % cut away corners of the hexagon
%                     p( p( :, 1 ).^2 + p( :, 2 ).^2 > ( diameter / ( 4 * cn ) )^2 * 4 / 3 + diameter^2 / 4, : ) = [];
%                 end
            elseif strcmp( rflag, 'square' )
                % find the closest square number to cnt
                rad = diameter / 2;
                %per = sqrt( pi * rad^2 / cnt ); % sqrt( area per ray )
                per = sqrt( diameter^2 / cnt ); % sqrt( area per ray )
                nr = ceil( rad / per ); % number of rays in each direction
                [ x, y ] = meshgrid( -nr * per : per : nr * per, -nr * per : per : nr * per );
                p( :, 1 ) = y( : );
                p( :, 2 ) = x( : );
                %p( p( :, 1 ).^2 + p( :, 2 ).^2 > rad^2, : ) = []; % remove corners
            elseif strcmp( rflag, 'pentile' )
                % generate a pentile grid
                dim = round( sqrt( cnt / 8 ) ); % number of cells in each dimension
                p = zeros( dim^2 * 8, 2 );
                p( 1, : ) = [ 0   0 ]; % green origin
                p( 2, : ) = [ .5   0 ]; % green right-bottom
                p( 3, : ) = [ 0   .5 ]; % green left-top
                p( 4, : ) = [ .5  .5 ]; % green right-top
                p( 5, : ) = [ .25 .25 ]; % red left-bottom
                p( 6, : ) = [ .75 .75 ]; % red right-top
                p( 7, : ) = [ .25 .75 ]; % blue left-top
                p( 8, : ) = [ .75 .25 ]; % blue right-bottom
                self.w( 1:4, 1 ) = 5300e-10;  % green wavelength of the OLED display
                self.color( 1:4, : ) = repmat( [ 0 1 0 ], 4, 1 );
                self.w( 5:6, 1 ) = 6200e-10;  % red wavelength of the OLED display
                self.color( 5:6, : ) = repmat( [ 1 0 0 ], 2, 1 );
                self.w( 7:8, 1 ) = 4580e-10;  % blue wavelength of the OLED display
                self.color( 7:8, : ) = repmat( [ 0 0 1 ], 2, 1 );
                
                for i = 0 : dim - 1
                    for j = 1 : dim
                        ind = 8 * ( i * dim + j );
                        p( ind + 1 : ind + 8, 1 ) = p( 1:8, 1 ) + j - ceil( dim/2 );
                        p( ind + 1 : ind + 8, 2 ) = p( 1:8, 2 ) + i - floor( dim/2 );
                        self.w( ind + 1 : ind + 8, 1 ) = self.w( 1:8 );
                        self.color( ind + 1 : ind + 8, : ) = self.color( 1:8, : );
                    end
                end
                p = p( 9 : end, : ); % remove the original cell
                self.w = self.w( 9 : end );
                self.color = self.color( 9 : end, : );
                p = p * diameter / 2 / max( p( :, 1 ) ); % scale the ray positions
                ind = p( :, 1 ).^2 + p( :, 2 ).^2 > diameter^2 / 4;
                p( ind, : ) = []; % leave rays only within the diameter
                self.w( ind ) = [];
                self.color( ind, : ) = [];
            elseif strcmp( rflag, 'sphere' ) % create a ray bundle with rays distributed in a spherical fashion around a source
                if ~strcmp( geometry, 'source' )
                    error( 'Sphere bundle pattern requires Source bundle geometry!' );
                end
                % find the closest square number to cnt
                n = round( sqrt( cnt ) );
                if rem( n, 2 ) == 1
                    n = n + 1;
                end
                [ x, y, z ] = sphere( n ); % create a spherical distribution of points along latitude and longtitude lines
                % rotate so that the first point in each const. z row faces the positive x direction and make the rows go along the y-change direction
                x = circshift( x, n/2 + 1, 2 )';
                y = circshift( y, n/2 + 1, 2 )'; 
                z = circshift( z, n/2 + 1, 2 )';
                self.n = [ x(:) y(:) z(:) ];
                self.cnt = size( self.n, 1 );
                self.r = repmat( pos, self.cnt, 1 );
                if norm( ax ) ~= 0
                    self.n = rodrigues_rot( self.n, ax, asin( norm( ax ) ) );
                end
            else
                error( [ 'Ray arrangement flag ' rflag ' is not defined!' ] );
            end
            
            if ~strcmp( rflag, 'sphere' )
                self.cnt = size( p, 1 );            
                p = [ zeros( self.cnt, 1 ) p ]; % add x-positions
                pos = repmat( pos, self.cnt, 1 );
                if norm( ax ) ~= 0
                    p = rodrigues_rot( p, ax, asin( norm( ax ) ) );
                end
                if strcmp( geometry, 'collimated' ) % parallel rays
                    % distribute over the area
                    self.r = pos + p;
                    dir = repmat( dir, self.cnt, 1 );
                    self.n = dir;
                elseif strcmp( geometry, 'source' ) || strcmp( geometry, 'source-Lambert' ) % assume p array at dir, source at pos.
                    self.r = pos;
                    self.n = p + repmat( dir, self.cnt, 1 );
                elseif strcmp( geometry, 'vergent' ) %
                    % distribute over the area
                    self.r = pos + p;
                    if diopter == 0 % the same as collimated
                        dir = repmat( dir, self.cnt, 1 );
                        self.n = dir;
                    else
                        self.n = p + repmat( 1000 * dir / diopter, self.cnt, 1 );
                    end
                    if diopter < 0
                        self.n = -self.n; % make ray normal point forward, as usual
                    end
                else
                    error( [ 'Source geometry' source ' is not defined!' ] );
                end
                % normalize directions
                self.n = self.n ./ repmat( sqrt( sum( self.n.^2, 2 ) ), 1, 3 );
            end
                
            if ~strcmp( rflag, 'pentile' )
                self.w = repmat( self.w, self.cnt, 1 );
                self.color = repmat( self.color, self.cnt, 1 );
            end
            self.nrefr = refrindx( self.w, glass );
            self.I = ones( self.cnt, 1 );
            if strcmp( geometry, 'source-Lambert' )
                self.I = self.I .* self.n( :, 1 ); % Lambertian source: I proportional to cos wrt source surface normal assumed to be [ 1 0 0 ]
            end
            
            self.att = ones( self.cnt, 1 );
            
            self.opl = zeros( self.cnt, 1 );
            self.gpl = zeros( self.cnt, 1 );
        end
        
            
        function draw( self, scale )
            if nargin == 0 || isempty( scale )
                scale = 1;
            end
            vis = self.I ~= 0;
            [ unique_colors, ~, ic ] = unique( self.color, 'rows' );
            nrms = scale * self.n;
            for i = 1 : size( unique_colors, 1 )
                cvis = vis & ( ic == i );
                quiver3( self.r( cvis, 1 ), self.r( cvis, 2 ), self.r( cvis, 3 ), ...
                         nrms( cvis, 1 ),   nrms( cvis, 2 ),   nrms( cvis, 3 ), ...
                         0, 'Color', unique_colors( i, : ), 'ShowArrowHead', 'off' );
            end
        end
         
        
        function [ rays_out, nrms ] = intersection( self, surf )
            % instantiate Rays object
            rays_out = self; % copy incoming rays
            
            switch class( surf )
                
                case { 'Aperture', 'Plane', 'Screen','ScreenGeneric' } % intersection with a plane
                    % distance to the plane along the ray
                    d = dot( repmat( surf.n, self.cnt, 1 ), repmat( surf.r, self.cnt, 1 ) - self.r, 2 ) ./ ...
                        dot( self.n, repmat( surf.n, self.cnt, 1 ), 2 ); 
                     
                    % calculate intersection vectors and normals
                    rinter = self.r + repmat( d, 1, 3 ) .* self.n;
                    nrms = repmat( surf.n, self.cnt, 1 );
                    

                    
                    
                    % bring surface to the default position
                    rtr = rinter - repmat( surf.r, self.cnt, 1 );
                    if surf.rotang ~= 0
                        rtr = rodrigues_rot( rtr, surf.rotax, -surf.rotang ); % rotate rays to the default plane orientation
                    end
                    
                    rays_out.r = rinter;
                    rinter = rtr;
                    if isa( surf, 'Screen' ) % calculate retinal image
                        wrong_dir = dot( nrms * sign( surf.R(1) ), self.n, 2 ) < 0;
                        self.I( wrong_dir ) = 0; % zero for the rays that point away from the screen for the image formation
                        rays_out.r( wrong_dir, : ) = Inf * rays_out.r( wrong_dir, : );
                        rays_out.opl( wrong_dir, : ) = NaN;                     
                        rays_out.gpl( wrong_dir, : ) = NaN; 
                        
                        surf.image = hist2( rtr( :, 2 ), rtr( :, 3 ), self.I, ...
                            linspace( -surf.w/2, surf.w/2, surf.wbins ), ...
                            linspace( -surf.h/2, surf.h/2, surf.hbins ) );
                        surf.image = flipud( surf.image ); % to get from matrix to image form
                        surf.image = fliplr( surf.image ); % because y-axis points to the left
                    end
                    
                    if isa( surf, 'ScreenGeneric' ) % store general info (opl,wf or tilt)
                        wrong_dir = dot( nrms * sign( surf.R(1) ), self.n, 2 ) < 0;
                        self.I( wrong_dir ) = 0; % zero for the rays that point away from the screen for the image formation
                        rays_out.r( wrong_dir, : ) = Inf * rays_out.r( wrong_dir, : );
                        rays_out.opl( wrong_dir, : ) = NaN;                     
                        rays_out.gpl( wrong_dir, : ) = NaN; 
                        
                        
                        switch surf.type
                            case 'wf'
                                z= 1/1000*(rays_out.opl-min(rays_out.opl))./self.w; 
                            case 'opl'
                                z= rays_out.opl;
                            case 'tilt'
                                z= acos(dot( nrms * sign( surf.R(1) ), self.n, 2 ));
                        end
                        
                        surf.raw = [rtr( :, 2 ), rtr( :, 3 ), z];
                        surf.image = hist2mean( rtr( :, 2 ), rtr( :, 3 ), z, ...
                            linspace( -surf.w/2, surf.w/2, surf.wbins ), ...
                            linspace( -surf.h/2, surf.h/2, surf.hbins ) );
                        surf.image = flipud( surf.image ); % to get from matrix to image form
                        surf.image = fliplr( surf.image ); % because y-axis points to the left
                    end
                    
                    
                    % handle rays that miss the element
                    out = [];
                    if isprop( surf, 'w' ) && ~isempty( surf.w ) && isprop( surf, 'h' ) && ~isempty( surf.h )
                        out =  rinter( :, 2 ) < -surf.w/2 | rinter( :, 2 ) > surf.w/2 | ...
                            rinter( :, 3 ) < -surf.h/2 | rinter( :, 3 ) > surf.h/2;
                    elseif isprop( surf, 'D' ) && ~isempty( surf.D )
                        if length( surf.D ) == 1
                            out = sum( rinter( :, 2:3 ).^2, 2 ) - 1e-12 > ( surf.D / 2 )^2;
                        elseif length( surf.D ) == 2
                            r2 = sum( rinter( :, 2:3 ).^2, 2 );
                            out = ( r2 + 1e-12 < ( surf.D(1) / 2 )^2 ) | ( r2 - 1e-12 > ( surf.D(2) / 2 )^2 );
                        elseif length( surf.D ) == 4
                            out =  rinter( :, 2 ) > -surf.D(1)/2 & rinter( :, 2 ) < surf.D(1)/2 & ...
                                   rinter( :, 3 ) > -surf.D(2)/2 & rinter( :, 3 ) < surf.D(2)/2 | ...
                                   rinter( :, 2 ) < -surf.D(3)/2 | rinter( :, 2 ) > surf.D(3)/2 | ...
                                   rinter( :, 3 ) < -surf.D(4)/2 | rinter( :, 3 ) > surf.D(4)/2;
                        end
                    end
                    rays_out.I( out ) = -1 * rays_out.I( out ); % mark for processing in the interaction function
                    
                    if isa( surf, 'Screen' )  % do not draw rays that missed the screen
                        rays_out.I( out ) = 0;
                        %rays_out.I( wrong_dir ) = 0; 
                        rays_out.r( out, : ) = Inf;
                    elseif isa( surf, 'Aperture' )
                        rays_out.I( ~out ) = 0; % block the rays
                    end
                    
                case { 'GeneralLens' 'AsphericLens' 'FresnelLens' 'ConeLens' 'CylinderLens' 'Lens' 'Sphere' 'Retina' } % intersection with a conical surface of rotation
                    % intersection between rays and the surface, also returns surface normals at the intersections
                    
                    % transform rays into the lens surface RF
                   
                    
                    r_in = self.r - repmat( surf.r, self.cnt, 1 ); % shift to RF with surface origin at [ 0 0 ]
                    
                    if surf.rotang ~= 0 % rotate so that the surface axis is along [1 0 0]
                        r_in = rodrigues_rot( r_in, surf.rotax, -surf.rotang ); % rotate rays to the default surface orientation
                        e = rodrigues_rot( self.n, surf.rotax, -surf.rotang );
                    else
                        e = self.n;
                    end
                    
                    if size( surf.R, 2 ) > 1 % asymmetric quadric, scale z-dimension to make the surface symmetric
                        sc = surf.R( 1 ) / surf.R( 2 );
                        r_in( :, 3 ) = r_in( :, 3 ) * sc;
                        e( :, 3 ) = e( :, 3 ) * sc;
                    end
                    
                    if isa( surf , 'GeneralLens' ) || isa( surf, 'AsphericLens' )
                        % minimize a measure of distance between a ray point and the surface
                        rinter = ones( self.cnt, 3 ); % init intersection vectors
                        outs = self.I == 0;
                        
                        
                        switch self.setget_solver_flag()
                            case {'default', 'precise'}
                              a_fun = @dist2; %avoid broadcasting in par loop
                              a_fval= 1e-8;
                            case 'precise_abs'
                              a_fun = @dist_abs; %avoid broadcasting in par loop
                              a_fval= 1e-4;
                            otherwise
                             warning('OptoMetrika:solverFlag',...
                                     'solverFlag must be ´default´ or ´precise´.' );   
                        end
                        
                        if exist( 'fminunc', 'file' ) % requires optimization toolbox
                            options = optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'Diagnostics', 'off');

                            if strcmp('precise',self.setget_solver_flag()) || ...
                                strcmp('precise_abs',self.setget_solver_flag()) 
                                options.StepTolerance = 1e-15;
                                options.OptimalityTolerance = 1e-15;
                                options.FiniteDifferenceStepSize = 100*eps;
                            end
                            
                            parfor i = 1 : self.cnt % run parallel computing
                                % for i = 1 : self.cnt % run normal computing
                                if outs( i ) == 0 % don't process lost rays
                                    [ d, fval ] = fminunc( a_fun, 20, options, r_in( i, : ), e( i, : ), surf );
                                    if fval > a_fval % didn't intersect with the surface
                                        outs( i ) = 1;
                                        rinter( i, : ) = Inf;
                                    else                                       
                                        rinter( i, : ) = r_in( i, : ) + e( i, : ) * d;
                                    end
                                end
                            end
                            rays_out.I( outs ) = 0;
                        else % no optimization toolbox
                            options = optimoptions( 'fminunc', 'MaxFunEvals', 2000, 'Display', 'off', 'Diagnostics', 'off' );
                            parfor i = 1 : self.cnt  % run parallel computing
                                if outs( i ) == 0 % don't process lost rays
                                    [ d, fval ] = fminsearch( a_fun, 0, options, r_in( i, : ), e( i, : ), surf );
                                    if fval > a_fval % didn't intersect with the surface
                                        outs( i ) = 1;
                                        rinter( i, : ) = Inf;
                                    else
                                        rinter( i, : ) = r_in( i, : ) + e( i, : ) * d;
                                        
                                    end
                                end
                            end
                        end % optimisation toolbox
                        
                        % get surface normals at the intersection points
                        en = surf.funch( rinter( :, 2 ), rinter( :, 3 ), surf.funca, 1 );
                        
                    elseif isa( surf , 'FresnelLens' ) % Fresnel lens
                        % find rays intersections with each Fresnel cone
                        rinter = Inf * ones( self.cnt, 3 ); % init intersection vectors
                        mem = [ rinter rinter ];
                        % rings = ones( self.cnt, 1 );        % ring indices of the rays
                        rings = zeros( self.cnt, 1 );        % ring indices of the rays
                        for i = 1 : surf.ncones
                            ren = [];
                            if i == 1
                                if length( surf.D ) == 1 || surf.D(1) == 0
                                    radin = 0;
                                else
                                    radin = 2 * surf.rad( 1, 1 ) - surf.rad( 2, 1 ); % take the inner radius assuming the same step
                                end
                            else
                                radin = surf.rad( i - 1, 1 );
                            end
                            if isempty( surf.vx ) || surf.R( i, 1 ) == 0 || isinf( surf.k( i ) ) % cone surface
                                [ in, rin ] = cone_intersection( r_in, e, radin, surf.rad( i, 2 ), surf.sag( i ), surf.the( i, 1 ), surf );
                            else % quadric surface
                                ring.k = surf.k( i );
                                ring.R = surf.R( i, 1 );
                                dims = size(r_in);
                                rin = repmat([ surf.vx( i ) + surf.sag( i ), 0, 0 ],dims(1),1) + ...
                                conic_intersection( r_in - repmat([ surf.vx( i ) + surf.sag( i ), 0, 0 ],dims(1),1), e, ring );
                                r2 = rin( :, 2 ).^2 + rin( :, 3 ).^2;
                                in = r2 >= radin.^2 & r2 < surf.rad( i, 2 ).^2;
                                if sum( in ) == 0
                                    continue;
                                end
                                % find normals
                                r2yz = r2 / ring.R^2; % distance to the lens center along the lens plane in units of lens R
                                if ring.k == -1 % parabola, special case
                                    c = 1 ./ sqrt( 1 + r2yz );
                                    s = sqrt( 1 - c.^2 );
                                else
                                    s = sqrt( r2yz ) ./ sqrt( 1 - ring.k * r2yz );
                                    c = sqrt( 1 - s.^2 );
                                end
                                s = -sign( ring.R ) * s; % sign of the transverse component to the ray determened by the lens curvature
                                th = atan2( rin( :, 3 ), rin( :, 2 ) ); % rotation angle to bring r into XZ plane
                                ren = [ c, s .* cos( th ), s .* sin( th ) ]; % make normal sign positive wrt ray
                                %figure, plot3( rin(:,1), rin(:,2), rin(:,3), 'b*' ), hold on, plot3( rin(in,1), rin(in,2), rin(in,3), 'r*' ), axis vis3d equal
                            end
                            d2old = sum( ( r_in - rinter ).^2, 2 ); % old distance between the ray start and the intersection, might be Inf
                            d2new = sum( ( r_in - rin ).^2, 2 ); % distance between the ray start and the intersection
                            closer = d2new < d2old;
                            vacant = ( rings == 0 ); % rays that didn't intersect any rings yet
                            if sum( vacant ) == 0
                                break;
                            end
                            in = in & closer & vacant;
                            mem( in, 1:3 ) = rin( in, : );
                            if ~isempty( ren )
                                mem( in, 4:6 ) = ren( in, : );
                            end
                            rays_out.I( in, : ) = self.I( in );
                            rings( in, 1 ) = i; % memorize ring indices for the intersecting rays
                        end
                        rinter = mem( :, 1:3 );
                        rings_hit = unique( rings );
 
                        % find normals
                        if isempty( surf.vx ) || surf.R( i, 1 ) == 0 || isinf( surf.k( i ) ) % cone surface
                            c = cos( surf.the( rings, 1 ) );
                            s = sin( surf.the( rings, 1 ) );
                            th = atan2( rinter( :, 3 ), rinter( :, 2 ) ); % rotation angle to bring r into XZ plane
                            en = repmat( sign( surf.the( rings, 1 ) ), 1, 3 ) .* [ s, -c .* cos( th ), -c .* sin( th ) ]; % make normal sign positive wrt ray
                        else
                            en = mem( :, 4:6 );
                        end
                        
                        % Correct for rays hitting at the very center of the lens, where the refraction is underfined. Assume that the refraction does not happen
                        central_rays = sum( rinter( :, 2:3 ).^2, 2 ) < 1e-20; % rays hitting the very center of the Fresnel lens
                        if sum( central_rays ) > 0
                            en( central_rays, : ) = self.n( central_rays, : );
                        end
                        
                        % find intersections with Fresnel walls
                        if surf.the( i, 2 ) == pi % vertical wall, cylinder
                            hits = false( self.cnt, 1 ); % record rays that hit the cylindrical walls here
                            mem = rinter;
                            for i = 1 : surf.ncones - 1
                                radin = surf.rad( i );
                                if i == 1
                                    if length( surf.D ) == 1
                                        rin = 0;
                                    else
                                        rin = 2 * surf.rad( 1, 1 ) - surf.rad( 2, 1 ); % take the inner radius assuming the same step
                                    end
                                    h = surf.sag( 1 ) + ( radin - rin ) / tan( surf.the( 1, 1 ) ) - surf.sag( 2 );
                                else
                                    h = surf.sag( i ) + ( radin - surf.rad( i - 1, 1 ) ) / tan( surf.the( i, 1 ) ) - surf.sag( i + 1 );
                                end
                                if radin > 0
                                    if h >= 0 % cylinder extends above the inner ring radius
                                        hsh = 0;
                                    else % cylinder extends below the inner ring radius
                                        h = -h ;
                                        hsh = h;
                                    end
                                    rsh = r_in - repmat( [ surf.sag( i + 1 ) - hsh, 0, 0 ], size( r_in, 1 ), 1 ); % shift rays origins back by sag
                                    [ in, rin ] = cylinder_intersection( rsh, e, 2 * radin, h, surf );
                                    %                                 if sum( in ) ~= 0
                                    %                                     figure, plot3( rinter( in, 1 ), rinter( in, 2 ), rinter( in, 3 ), 'r*' ), hold on, plot3( rsh( in, 1 ), rsh( in, 2 ), rsh( in, 3 ), 'ko' ), plot3( rin( in, 1 ), rin( in, 2 ), rin( in, 3 ), 'b*' ), axis equal vis3d
                                    %                                 end
                                    d2old = sum( ( r_in - rinter ).^2, 2 ); % old distance between the ray start and the intersection, might be Inf
                                    d2new = sum( ( rsh - rin ).^2, 2 ); % distance between the ray start and the intersection
                                    closer = d2new < d2old;
                                    in = in & closer; % only consider rays that missed Fresnel cones and are closer than previous intersections
                                    mem( in, : ) = rin( in, : ) + repmat( [ surf.sag( i + 1 ) - hsh, 0, 0 ], sum( in ), 1 );
                                    if surf.walls == 1 % wall from the same material, consider rays coming through
                                        rays_out.I( in, : ) = self.I( in );
                                    else % walls covered with soot
                                        rays_out.I( in, : ) = 0;
                                    end
                                    hits( in ) = 1;
                                end
                            end
                            rinter = mem;
                            % find normals
                            if sum( hits ) > 0
                                th = atan2( rinter( hits, 3 ), rinter( hits, 2 ) ); % rotation angle to bring r into XZ plane
                                en( hits, : ) = [ zeros( size( th ) ), cos( th ), sin( th ) ];
                            end
                        else % non-vertical wall, cone
                            mem = rinter;
                            for i = 1 : surf.ncones - 1
                                if i == 1
                                    if length( surf.D ) == 1 || surf.D(1) == 0
                                        radin = 0;
                                    else
                                        radin = 2 * surf.rad( 1, 1 ) - surf.rad( 2, 1 ); % take the inner radius assuming the same step
                                    end
                                else
                                    radin = surf.rad( i - 1, 1 );
                                end
                                if surf.rad( i, 2 ) < surf.rad( i, 1 )
                                    [ in, rin ] = cone_intersection( r_in, e, surf.rad( i, 2 ), surf.rad( i, 1 ), ...
                                        surf.sag( i ) + ( surf.rad( i, 2 ) - radin ) ./ tan( surf.the( i, 1 ) ), surf.the( i, 2 ), surf );
                                else
                                    [ in, rin ] = cone_intersection( r_in, e, surf.rad( i, 1 ), surf.rad( i, 2 ), ...
                                        surf.sag( i + 1 ), surf.the( i, 2 ), surf );
                                end
                                if ~isreal( rin )
                                    rin;
                                end
                                
                                d2old = sum( ( r_in - rinter ).^2, 2 ); % old distance between the ray start and the intersection, might be Inf
                                d2new = sum( ( r_in - rin ).^2, 2 ); % distance between the ray start and the intersection
                                closer = d2new < d2old;
                                in = in & closer;
                                if sum( in ) > 0
                                    mem( in, : ) = rin( in, : );
                                    if surf.walls == 1 % wall from the same material, consider rays coming through
                                        rays_out.I( in, : ) = self.I( in );
                                    else % walls covered with soot
                                        rays_out.I( in, : ) = 0;
                                    end
                                    % find normals
                                    c = cos( surf.the( i, 2 ) );
                                    s = sin( surf.the( i, 2 ) );
                                    th = atan2( rin( in, 3 ), rin( in, 2 ) ); % rotation angle to bring r into XZ plane
                                    en( in, : ) = repmat( sign( surf.the( i, 2 ) ), length( th ), 3 ) .* [ repmat( s, length( th ), 1 ), -c .* cos( th ), -c .* sin( th ) ]; % make normal sign positive wrt ray
                                end
                            end
                            rinter = mem;
                        end
                        
%                         figure, plot3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ), '*' ), hold on, ...
%                         quiver3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ), en( :, 1 ), en( :, 2 ), en( :, 3 ), 5 ), axis equal vis3d;
                      
                    elseif isa( surf, 'ConeLens' ) % cone lens
                        % find intersections
                        rinter = Inf * ones( self.cnt, 3 ); % init intersection vectors
                        [ in, rin ] = cone_intersection( r_in, e, surf.rad(1), surf.rad(2), 0, surf.the, surf );
                        rinter( in, : ) = rin( in, : );
                        rays_out.I( in, : ) = self.I( in );
                        % find normals
                        c = cos( surf.the );
                        s = sin( surf.the );
                        th = atan2( rinter( :, 3 ), rinter( :, 2 ) ); % rotation angle to bring r into XZ plane
                        en = sign( surf.the ) * [ s * ones( size( th ) ), -c * cos( th ), -c * sin( th ) ]; % make normal sign positive wrt ray
                                                
                    elseif isa( surf, 'CylinderLens' )
                        % find intersection
                        rinter = Inf * ones( self.cnt, 3 ); % init intersection vectors
                        [ in, rin ] = cylinder_intersection( r_in, e, surf.D(1), surf.h, surf );
                        rinter( in, : ) = rin( in, : );
                        rays_out.I( in, : ) = self.I( in );
                        % find normals
                        th = atan2( rinter( :, 3 ), rinter( :, 2 ) ); % rotation angle to bring r into XZ plane
                        en = [ zeros( size( th ) ), cos( th ), sin( th ) ];
       
                    elseif isa( surf, 'Sphere' )   
                        [rinter] = sphere_intersection( r_in, e, surf );  
                        

                        en =  rinter ;  
                        % normal vector coincident with
                         % intersection position
                        
                    else % conic lens
                        [rinter] = conic_intersection( r_in, e, surf );                        
                        
                        % find normals
                        r2yz = ( rinter( :, 2 ).^2 + rinter( :, 3 ).^2 ) / surf.R(1)^2; % distance to the lens center along the lens plane in units of lens R
                        if surf.k == -1 % parabola, special case
                            c = 1 ./ sqrt( 1 + r2yz );
                            s = sqrt( 1 - c.^2 );
                        else
                            s = sqrt( r2yz ) ./ sqrt( 1 - surf.k * r2yz );
                            c = sqrt( 1 - s.^2 );
                        end
                        %if R > 0 % add corrugations
                        %amp = 0.001;
                        %per = 1;
                        %c = 1 ./ sqrt( 1 + ( sqrt( r2yz ) ./ sqrt( 1 - ( 1 + k ) * r2yz ) + 2 * pi * amp / per * sin( 2 * pi / per * R * sqrt( r2yz ) ) ).^2 );
                        %s = sqrt( 1 - c.^2 );
                        %end
                        s = -sign( surf.R(1) ) * s; % sign of the transverse component to the ray determined by the lens curvature
                        th = atan2( rinter( :, 3 ), rinter( :, 2 ) ); % rotation angle to bring r into XZ plane
                        en = [ c, s .* cos( th ), s .* sin( th ) ]; % make normal sign positive wrt ray
                        
                        if isa( surf, 'Retina' ) % calculate retinal image
                            wrong_dir = dot( en, e, 2 ) < 0;
                            self.I( wrong_dir ) = 0; % zero for the rays that point away from the screen for the image formation
                            rinter( wrong_dir, : ) = NaN;
                            rinter( self.I == 0, : ) = NaN;
                            rinter( sqrt( sum( rinter.^2, 2 ) ) > realmax / 2, : ) = NaN;
                            % scale to a unit spherical surphace
                            rtr = rinter .* ( 1 + surf.k ) / -surf.R(1);
                            rtr( :, 1 ) = rtr( :, 1 ) + 1;
                            [ az, el ] = cart2sph( rtr( :, 1 ), rtr( :, 2 ), rtr( :, 3 ) ); % YZX to account for Optometrika's coordinate system
                            maz = max( abs( az ) );
                            mel = max( abs( el ) );
                            md = max( maz, mel );
                            if isfinite( md )
                                surf.image = hist2( az, el, self.I, ...
                                    linspace( -md, md, surf.azbins ), ...
                                    linspace( -md, md, surf.elbins ) );
                            end
                        end %is a retina
                    end %is general lens, aspheric lens,..
                                        
                    % handle rays that miss the element
                    out = [];
                    if isa( surf, 'ConeLens' ) || isa( surf, 'CylinderLens' )
                        out = ~in;
                    else
                        if isprop( surf, 'w' ) && ~isempty( surf.w ) && isprop( surf, 'h' ) && ~isempty( surf.h )
                            out =  rinter( :, 2 ) < -surf.w/2 | rinter( :, 2 ) > surf.w/2 | ...
                                rinter( :, 3 ) < -surf.h/2 | rinter( :, 3 ) > surf.h/2;
                        elseif isprop( surf, 'D' ) && ~isempty( surf.D )
                            if length( surf.D ) == 1
                                out = sum( rinter( :, 2:3 ).^2, 2 ) - 1e-12 > ( surf.D / 2 )^2;
                            else
                                r2 = sum( rinter( :, 2:3 ).^2, 2 );
                                out = isnan( r2 ) | ( r2 + 1e-12 < ( surf.D(1) / 2 )^2 ) | ( r2 - 1e-12 > ( surf.D(2) / 2 )^2 );
                            end
                        end
                    end
                    
                    if isa( surf, 'Retina' )  % do not draw rays that missed the screen
                        rays_out.I( out ) = 0;
                        rays_out.r( out, : ) = Inf;
                    else
                        rays_out.I( out ) = -1 * rays_out.I( out ); % mark for processing in the interaction function
                    end
                    
                    % return to the original RF
                    if size( surf.R, 2 ) > 1 % asymmetric quadric, unscale the z-dimension
                        sc = surf.R( 1 ) / surf.R( 2 );
                        rinter( :, 3 ) = rinter( :, 3 ) / sc;
                        en( :, 3 ) = en( :, 3 ) * sc; % normals transform as one-forms rather than vectors. Hence, divide by the scaling factor
                    end                  
                    if surf.rotang ~= 0 % needs rotation
                        rays_out.r = rodrigues_rot( rinter, surf.rotax, surf.rotang );
                        nrms = rodrigues_rot( en, surf.rotax, surf.rotang );
                    else
                        rays_out.r = rinter;
                        nrms = en;
                    end
                    nrms = nrms ./ repmat( sqrt( sum( nrms.^2, 2 ) ), 1, 3 );
%                     if ~isreal( en )
%                         % rays_out.I( imag( sum( en, 2 ) ) ~= 0 ) = 0;
%                         figure, plot3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ), '*' ), hold on,
%                         quiver3( rinter( :, 1 ), rinter( :, 2 ), rinter( :, 3 ), en( :, 1 ), en( :, 2 ), en( :, 3 ), 5 ),
%                         axis equal vis3d; xlabel( 'x' ), ylabel( 'y' ), zlabel( 'z' );
%                     end
                    rays_out.r = rays_out.r + repmat( surf.r, self.cnt, 1 );
                    
                otherwise
                    error( [ 'Surface ' class( surf ) ' is not defined!' ] );
            end
            
        end
         
        
        function rays_out = interaction( self, surf, out_fl )
            % INTERACTION calculates rays properties after interacting with
            % a Surface
            
            % find intersections and set outcoming rays starting points
            [ rays_out, nrms ] = self.intersection( surf );
                    
        
           %gpl1= dot( rays_out.r-self.r, self.n, 2 );
          rays_out.gpl =  dot( rays_out.r-self.r, self.n, 2 );
          rays_out.opl =  self.opl+rays_out.gpl.*self.nrefr;
            
            
            miss = rays_out.I < 0; % indices of the rays

            rays_out.gpl(miss)=NaN;

            med1 = surf.glass{1};
            med2 = surf.glass{2};
            
            
            % determine refractive indices before and after the surface 
            cs1 = dot( nrms, self.n, 2 ); % cosine between the ray direction and the surface direction
            opp_rays = cs1 < 0; %self.nrefr == refrindx( self.w, med2 ); %cs1 < 0; % rays hitting the surface from the opposite direction          
            old_refr( ~opp_rays ) = refrindx( self.w( ~opp_rays ), med1 ); % refractive index before the surface
            old_refr(  opp_rays ) = refrindx( self.w(  opp_rays ), med2 ); % refractive index before the surface
            if strcmp( med2, 'mirror' )
                new_refr = self.nrefr; % refractive index after the surface
                %old_refr = new_refr;
            elseif  strcmp( med1, 'mirror' )
                new_refr = refrindx( self.w, med2 ); % refractive index after the surface
                old_refr = new_refr;
            else
                new_refr( ~opp_rays ) = refrindx( self.w( ~opp_rays ), med2 ); % refractive index after the surface
                new_refr(  opp_rays ) = refrindx( self.w(  opp_rays ), med1 ); % refractive index after the surface
            end
            
            if size( new_refr, 1 ) < size( new_refr, 2 )
                new_refr = new_refr';
            end
            if size( old_refr, 1 ) < size( old_refr, 2 )
                old_refr = old_refr';
            end
            
              % Check the refractive index is correct
            if ~(self.nrefr == old_refr ) & ...
                ~isa( surf, 'Retina' ) & ~isa( surf, 'Screen' ) & ...
                ~isa( surf, 'ScreenGeneric' ) & ~isa( surf, 'Aperture' ) 
            
				warning('OptoMetrika:refrIndexMismatch',...
                    ['Refractive Index of Ray and Optical element mismatch. \n' ...
                     'Was changed to Value expected at element entry.\n' ]);
                self.nrefr = old_refr; 
                rays_out.nrefr = new_refr;
            end

            % calculate refraction
            switch( class( surf ) )
                case { 'Aperture', 'Screen', 'Retina','ScreenGeneric' } 
                case { 'GeneralLens' 'AsphericLens' 'FresnelLens' 'ConeLens' 'CylinderLens' 'Plane' 'Sphere' 'Lens' }
                    % calculate refraction (Snell's law)

                    inside_already = ( ~miss ) & ( abs( rays_out.nrefr - old_refr ) > 1e-12 ); % rays that are already inside the surface (entered it previously)
                    rays_out.nrefr( ~miss ) = new_refr( ~miss ); % change refractive index of the rays that crossed the surface
                    if sum( inside_already ) ~= 0 % second intersections in a cylinder
                        rays_out.nrefr( inside_already ) = old_refr( inside_already ); % use old refractive index for those rays that are crossing the second surface
                    end
            
%                     if isa( surf, 'FresnelLens' ) % remove rays coming to the Fresnel surface from the inside (through the cylindrical walls).
%                         bads = cs1 < 0;
%                         rays_out.I( bads ) = 0;
%                         rays_out.r( bads, : ) = Inf * rays_out.r( bads, : );
%                     end
                    
                    if strcmp( med1, 'mirror' ) || strcmp( med2, 'mirror' ) % if a mirror
                        rays_out.n = self.n - 2 * repmat( cs1, 1, 3 ) .* nrms; % Snell's law of reflection
                        %rays_out.nrefr = refrindx( self.w, med1 ); % refractive index before the surface
                        if strcmp( med1, 'mirror' ) && strcmp( med2, 'air' ) % mirror facing away
                            rays_out.I( cs1 > 0 & ~miss ) = 0; % zero rays hitting such mirror from the back
                        elseif strcmp( med1, 'air' ) && strcmp( med2, 'mirror' ) % mirror facing toward me
                            rays_out.I( cs1 < 0 & ~miss ) = 0; % zero rays hitting such mirror from the back
                        end
                    elseif strcmp( med1, 'soot' ) || strcmp( med2, 'soot' ) % opaque black
                        rays_out.I( ~miss ) = 0; % zero rays that hit the element
                    else % transparent surface
                        rn = self.nrefr ./ rays_out.nrefr; % ratio of in and out refractive indices
                        cs2 = sqrt( 1 - rn.^2 .* ( 1 - cs1.^2 ) );
                        rays_out.n = repmat( rn, 1, 3 ) .* self.n - repmat( rn .* cs1 - sign( cs1 ) .* cs2, 1, 3 ) .* nrms; % refracted direction
                        tmp = cs1;
                        cs1( opp_rays ) = -cs1( opp_rays );
                        
                        % calculate transmitted intensity (Fresnel formulas)
                        rs = ( rn .* cs1 - cs2 ) ./ ( rn .* cs1 + cs2 );
                        rp = ( cs1 - rn .* cs2 ) ./ ( cs1 + rn .* cs2 );
                        refraction_loss = ( abs( rs ).^2 + abs( rp ).^2 ) / 2;
                        % handle total internal reflection
                        tot = find(imag( cs2 ) ~= 0);
                        % rays_out.n( tot, : ) = 0; % zero direction for such rays
                        if(~isempty(tot))   
                            rays_out.n( tot, : ) = self.n( tot, : ) - 2 * repmat( tmp( tot ), 1, 3 ) .* nrms( tot, : ); % Snell's law of reflection
                            refraction_loss( tot ) = 0; %1;
                            rays_out.nrefr( tot ) = refrindx( self.w( tot ), med1 ); % refractive index before the surface
                        end
                        
                        rays_out.I( ~miss ) = ( 1 - refraction_loss( ~miss ) ) .* rays_out.I( ~miss ); % intensity of the outcoming rays
                    end                     
                otherwise
                    error( [ 'Surface ' class( surf ) ' is not defined!' ] );
            end
            
            % process rays that missed the element
            if out_fl == 0 || strcmp( med1, 'soot' ) || strcmp( med2, 'soot' ) % if tracing rays missing elements or for apertures
                % use the original rays here
                rays_out.I( miss ) = self.I( miss );
                rays_out.r( miss, : ) = self.r( miss, : );
                rays_out.n( miss, : ) = self.n( miss, : );
            else
                % default, exclude such rays
                rays_out.I( miss ) = 0;
                rays_out.r( miss, : ) = Inf;
            end
            rays_out.I( isnan( rays_out.I ) ) = 0;
            %rays_out.I( rays_out.n( :, 1 ) < 0 ) = 0; % zero rays that point back to the source
        end
        
        
        function rc = copy( self )
            rc = Rays;  % initialize the class instance
            rc.r = self.r; % a matrix of ray starting positions
            rc.n = self.n; % a matrix of ray directions
            rc.w = self.w;  % a vector of ray wavelengths
            rc.I = self.I;         % a vector of ray intensities
            rc.nrefr = self.nrefr;    % a vector of current refractive indices
            rc.att = self.att;       % a vector of ray attenuations
            rc.color = self.color;   % color to draw the bundle rays
            rc.cnt = self.cnt;      % number of rays in the bundle   
            rc.opl = self.opl;      
            rc.gpl = self.gpl;     
        end
        
        
        function self = append( self, rays )
            % append rays to the current bundle
            self.r = [ self.r; rays.r ];
            self.n = [ self.n; rays.n ];
            self.w = [ self.w; rays.w ];
            self.I = [ self.I; rays.I ];
            self.opl = [ self.opl; rays.opl ];
            self.gpl = [ self.gpl; rays.gpl ];
            self.nrefr = [ self.nrefr; rays.nrefr ];
            self.att = [ self.att; rays.att ];
            self.color = [ self.color; rays.color ];
            self.cnt = self.cnt + rays.cnt;                
        end
        
        function rays = subset( self, inds )
            % pick a subset of rays defined by inds in the current bundle
            rays = Rays; % allocate an instance of rays
            rays.r = self.r( inds, : );
            rays.n = self.n( inds, : );
            rays.w = self.w( inds, : );
            rays.I = self.I( inds, : );
            rays.gpl = self.gpl( inds, : );
            rays.opl = self.opl( inds, : );
            rays.nrefr = self.nrefr( inds, : );
            rays.att = self.att( inds, : );
            rays.color = self.color( inds, : );
            rays.cnt = length( inds );
        end
        
        function self = truncate( self )
            % remove rays with zero intensity
            ind = self.I == 0;
            self.r( ind, : ) = [];
            self.n( ind, : ) = [];
            self.w( ind, : ) = [];
            self.I( ind, : ) = [];
            self.opl( ind, : ) = [];
            self.gpl( ind, : ) = [];
            self.color( ind, : ) = [];
            self.nrefr( ind, : ) = [];
            self.att( ind, : ) = [];
            self.cnt = self.cnt - sum( ind );
        end
                            
        function [ av, dv, nrays ] = stat( self )
            % calculate mean and standard deviation of the rays startingpoints
            vis = self.I ~= 0; % visible rays
            norm = sum( self.I( vis ) );
            av = sum( repmat( self.I( vis ), 1, 3 ) .* self.r( vis, : ) ) ./ norm;
            dv = sqrt( sum( self.I( vis ) .* sum( ( self.r( vis, : ) - repmat( av, sum( vis ), 1 ) ).^2, 2 ) ) ./ norm );
            nrays = sum( vis );
        end

        function [ x0, cv, ax, ang, nrays ] = stat_ellipse( self )
            % calculate parameters of an ellips covering the projection of
            % the rays startingpoints onto YZ plane
            vis = self.I ~= 0; % visible rays
            nrays = sum( vis );
            [ x0, cv, ax, ang ] = ellipse_fit( self.r( vis, 2:3 ) );
            %[ x0, cv, ax, ang ] = ellipse_fit( self.r( vis, : ) );
        end

        function [ mu, sigma, lambda, angle, nrays ] = stat_exGaussian( self )
            % calculate parameters of the exponentially modified Gaussian distribution (EMG) fitting the rays startingpoints
            vis = self.I ~= 0; % visible rays
            nrays = sum( vis );
            [ mu, sigma, lambda, angle ] = exGaussian_fit2D( self.r( vis, 2:3 ) );
        end
        
        function [ av, dv ] = stat_sph( self, sph_pos )
            % calculate mean and standard deviation of the rays
            % startingpoints in spherical coordinates (e.g., on a retina)
            % returns average and std for [ azimuth, elevation, radius ]
            rs = self.r - sph_pos; % coordinates wrt sphere center
            sph = cart2sph( rs( :, 1 ), rs( :, 2 ), rs( :, 3 ) ); % [ az el r ]
            vis = self.I ~= 0; % visible rays
            norm = sum( self.I( vis ) );
            av = sum( repmat( self.I( vis ), 1, 3 ) .* sph( vis, 1:2 ) ) ./ norm;
            dv = sqrt( sum( self.I( vis ) .* sum( ( sph( vis, 1:2 ) - repmat( av, sum( vis ), 1 ) ).^2, 2 ) ) ./ norm );
        end
        
        function r2 = dist2rays( self, p )
            t = self.r - repmat( p, self.cnt, 1 ); % vectors from the point p to the rays origins
            proj = t * self.n'; % projections of the vectors t onto the ray
            r2 = sum( ( t - repmat( proj, 1, 3 ) .* self.n ).^2 );
        end
        
         function [ ndev, ndev_max ] = stat_div( self )
            % calculate rays divergence, avg and max
            
            ind = (self.I ~= 0);
            sn = self.n( ind, : );
            si = self.I( ind, : );
            repI = repmat( si , 1, 3 );
            nav = sum( sn .* repI, 1 ); % average bundle direction
            nav = nav / sqrt( sum( nav.^2, 2 ) ); % normalize the average direction vector
            nangle= acos(dot(self.n,repmat( nav, self.cnt,1),2 ));
            
            ndev= mean(abs(nangle(ind)));
            ndev_max= max(abs(nangle(ind)));
        end
        
        
        function [ f, ff ] = focal_point( self, flag )
            if nargin < 2
                flag = 0;
            end
            ind = self.I ~= 0;
            sn = self.n( ind, : );
            sr = self.r( ind, : );
            si = self.I( ind, : );
            repI = repmat( si , 1, 3 );
            nav = sum( sn .* repI, 1 ); % average bundle direction
            nav = nav / sqrt( sum( nav.^2, 2 ) ); % normalize the average direction vector
            tmp = repmat( nav, size( sr, 1 ), 1 );
            osr = tmp .* repmat( dot( sr, tmp, 2 ), 1, 3 );
            sr = sr - osr; % leave only the r component orthogonal to the average bundle direction.
            osr = sum( osr .* repI ) / sum( self.I ); % bundle origin along the bundle direction
            rav = sum( sr .* repI ) / sum( self.I ); % average bundle origin 

            if flag == 2 % search for the plane with the smallest cross-section
                scnt = sum( ind );
                if exist( 'fminunc', 'file' ) % requires optimization toolbox
                    options = optimoptions( 'fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'Diagnostics', 'off' );
                    plr = rav + nav * fminunc( @scatterRayPlaneIntersect, 0, options, sr, sn, si, rav, nav ); % optimal plane position
                else % no optimization toolbox
                    options = optimoptions( 'fminunc', 'MaxFunEvals', 2000, 'Display', 'off', 'Diagnostics', 'off' );
                    plr = rav + nav * fminsearch( @scatterRayPlaneIntersect, 0, options, sr, sn, si, rav, nav );
                end
                d = dot( repmat( nav, scnt, 1 ), repmat( plr, scnt, 1 ) - sr, 2 ) ./ ...
                    dot( sn, repmat( nav, scnt, 1 ), 2 );
                % calculate intersection vectors
                rinter = sr + repmat( d, 1, 3 ) .* sn; % intersection vectors
                f = mean( rinter ); % assign focus to the mean of the intersection points
                ff = mean( sqrt( sum( ( rinter - repmat( f, size( rinter, 1 ), 1 ) ).^2, 2 ) ) );
            else % use average cosine
                dr = sr - repmat( rav, size( sr, 1 ), 1 );
                ndr = sqrt( sum( dr.^2, 2 ) );
                drn = dr ./ repmat( ndr, 1, 3 ); % normalize the difference of origins vector
                dp = dot( sn, drn, 2 );
                if flag == 1  % only consider rays in the bundle which are diverging, the ones an eye can focus
                    ind = sign( dp ) > 0;
                    if sum( ind ) == 0
                        ind = ~ind; % consider the other half then
                    end
                else
                    ind = ones( size( si ) );
                end
                ssi = sum( si( ind ) ); % normalization factor for the weighted average
                cs = -sum( dp( ind ) .* si( ind ) ) / ssi; % mean cosine of the convergence angles
                d = sum( ndr( ind ) .* si( ind ) ) / ssi * sqrt( 1 - cs^2 ) / cs; % distance to the focus along the mean direction
                f = osr + d * nav; % focal point
                % calculate ray positions at the focal plane
                rf = dr + repmat( d ./ sqrt( 1 - dp.^2 ), 1, 3 ) .* sn;
                arf = mean( rf ); % average vector
                ff = sum( si .* sqrt( sum( ( rf - repmat( arf, size( rf, 1 ), 1 ) ).^2, 2 ) ) ) / ssi; % weighted average bundle size on the focal plane
            end               
        end
        
        function [ mtf, fr ] = MTF( self, dist, fr )
            if nargin < 3
                fr = linspace( 0, 50, 100 ); % frequencies to calculate MTF at
            end           
            ind = self.I ~= 0;
            sn = self.n( ind, : );
            sr = self.r( ind, : );
            si = self.I( ind, : );
            repI = repmat( si , 1, 3 );
            nav = sum( sn .* repI, 1 ); % average bundle direction
            nav = nav / sqrt( sum( nav.^2, 2 ) ); % normalize the average direction vector
            tmp = repmat( nav, size( sr, 1 ), 1 );
            osr = tmp .* repmat( dot( sr, tmp, 2 ), 1, 3 );
            sr = sr - osr; % leave only the r component orthogonal to the average bundle direction.
            osr = sum( osr .* repI ) / sum( self.I ); % bundle origin along the bundle direction
            rav = sum( sr .* repI ) / sum( self.I ); % average bundle origin
            
            [ ~, rinter ] = scatterRayPlaneIntersect( dist, sr, sn, si, rav, nav );
            % rotate the intersection vectors for the average to face in the x-direction
            rv = cross( nav, [ 1 0 0 ] );
            ra = asin( norm( rv ) );
            rinrot = rodrigues_rot( rinter, rv, ra );
            yz = rinrot( :, 2:3 );
            lsfY = 180 / pi * atan( ( yz( :, 1 ) - mean( yz( :, 1 ) ) ) / dist ); % convert to degrees
            lsfZ = 180 / pi * atan( ( yz( :, 2 ) - mean( yz( :, 2 ) ) ) / dist ); % convert to degrees
            %figure, scatter( lsfY, lsfZ, '*' ), axis equal;
            mtf = sum( cos( 2 * pi * lsfZ * fr ) );
            mtf = mtf / mtf( 1 );
        end
   end
end


