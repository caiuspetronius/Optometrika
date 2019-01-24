classdef Eye < Bench
    % EYE implements human eye optics
    %   The human eye model is described in Escudero-Sanz & Navarro, 
    %   "Off-axis aberrations of a wide-angle schematic eye model", 
    %   JOSA A, 16(8), 1881-1891 (1999). Also, Dubbelman M, Van der 
    %   Heijde GL. "The shape of the aging human lens: curvature, 
    %   equivalent refractive index and the lens paradox." Vision Res. 
    %   2001;41:1867?1877 and A. V. Goncharov and C. Dainty. Wide-field 
    %   schematic eye models with gradient-index lens. J Opt Soc Am A 
    %   Opt Image Sci Vis, 24(8):2157?74, 2007. were used.
    %   
    %   Following  G. K. Von Noorden and E. C. Campos. "Binocular 
    %   vision and ocular motility: theory and management of strabismus." 
    %   Gunter K. 6th ed. St. Louis: CV Mosby, 2002, the eye center of 
    %   rotation was taken to be 13.3 mm behind the corneal apex, which 
    %   puts it 1.34 mm behind the center of the eye (total eye depth is
    %   23.93 here). High precision is immaterial since eye center
    %   of rotation actually moves by as much as 1 mm as the eye rotates,
    %   so the idea of the eye center of rotation is just an approximation.
    %   The retina shape was taken to be an oblate spheroid according to 
    %   D. A. Atchison et al. "Shape of the retinal surface in emmetropia 
    %   and myopia". Invest Ophthalmol Vis Sci, 46(8):2698?707, Aug 2005.
    %   The retinal spheroid and lens slants were ignored here.
    %
    %   Lens accomodation is modeled by its diameter variation assuming 
    %   that the lens volume is constant for all diameters. This assumption
    %   was based on Hermans et al., "Constant Volume of the Human Lens and 
    %   Decrease in Surface Area of the Capsular Bag during Accommodation: 
    %   An MRI and Scheimp?ug Study", Investigative Ophthalmology & Visual 
    %   Science 50(1), 281-289 (2009).
    %
    %   The back lens surface was modeled by a paraboloid of revolution 
    %   x = 1/(2 R) ( y^2 + z^2 ), its (constant) volume V is 
    %   given by pi/8 D^2 h, where h is the paraboloid height. Hence, 
    %   h(D) = 8 V / (pi D^2) and R(D) = pi D^4 / (64 V). The front lens 
    %   surface was modeled by a hyperboloid of revolution given by
    %   x = R/(1+k) (1 - sqrt( 1 - a( y^2 + z^2)/R^2 ) ). The volume 
    %   is pi/8 D^2 h (1 - h/(6 R/(1+k) + 3h) ), where h is the height of
    %   of the hyperboloid from it apex to the cut of diameter D. 
    %   Corresponding h(D) and R(D) formulas were obtained by solving a cubic 
    %   equation, and its closed-from solution is implemented here.
    %
    % Member functions:
    %
    % e = Eye( lens_diameter, pupil_diameter )  - constructor function
    % INPUT:
    %   lens_diameter - eye lens diameter, mm, specifies the eye's
    %   accomodative state, defaults to no accommodation
    %   pupil_diameter - eye pupil diameter, mm, defaults to 3
    % OUTPUT:
    %   e - eye object
    %
    % e = e.Lens( lens_diameter ) - modify the eye's accommodative state
    % INPUT:
    %   lens_diameter = new lens diameter, mm
    % OUTPUT:
    %   e - eye object with the lens diameter modified
    %
    % V = e.eye_lens_vol( lens_diameter ) - returns the eye lens's volume
    % INPUT:
    %   lens_diameter - eye lens diameter, mm
    % OUTPUT:
    %   V - volume of the lens, mm^3
    %
    % e.rotate( rot_axis, rot_angle ) - rotate the eye e with all its elements
    % INPUT:
    %   rot_axis - 1x3 vector defining the rotation axis
    %   rot_angle - rotation angle (radians)
    %   
    % e.translate( tr_vec ) - translate the eye e with all its elements
    % INPUT:
    %   tr_vec - 1x3 translation vector
    %
    % e.display() - displays the eye e information
    %
    % e.draw() - draws the eye e in the current axes
    %
    % Copyright: Yury Petrov, 2016
    %
    
    properties ( Constant )
        % All distances wrt the eye center of rotation
        % cornea
        Cornea1D = 11     % front surface diameter, mm
        Cornea1R = 7.76   % front surface radius of curvature, mm (30-year old, Goncharov et al. 2009)
        Cornea1k = -0.10  % front conic constant (30-year old, Goncharov et al. 2009, Escudero-Sanz et al. 1999 uses -0.26
        Cornea1x = -13.3  % surface position wrt eye center of rotation
        Cornea2D = 10.5     % back surface diameter
        Cornea2R = 6.52   % back surface radius of curvature (30-year old, Goncharov et al. 2009), 6.50 (Escudero-Sanz et al. 1999)
        Cornea2k = -0.30  % back conic constant (30-year old, Goncharov et al. 2009), 0 (Escudero-Sanz et al. 1999) 
        Cornea2x = -12.75 % back surface position, cornea thickness is .55 (both Goncharov and Escudero-Sanz)
        % pupil
        PupilDout = 11    % outer pupil diameter 
        Pupilx   = -10.25 % -9.85  % Aqueous thickness is 3.05 (Goncharov, 2009 30-year old and also Escudero-Sanz & Navarro). 
                          % I shifted pupil plane .15 mm anterior to
                          % prevent front lens surface getting into the
                          % pupil surface for the smallest pupil (2mm) and
                          % the strongest accommodation (7D)
        % lens
        Lens1k   = -3.1316 % strongly hyperbolic
        Lens2k   = -1.0   % parabolic surface
        LensMx   = -8.55  % x position of the lens hyp|par boundary (calculated from Goncharov et al, 2009)
        LensV1   = 53.49  % constant volume of the lens hyperbolic part
        LensV2   = 106.66 % constant volume of the lens parabolic part
        % retina
        RetinaR  = -12.94  % radius of curvature at the apex or 1/2 retina length along optic axis
        Retinak  = 0.275  % emmetropic retina is oblong
        Retinax  = 10.6203 % retina apex position, this puts equatorial plane at .47
        % RetinaRotz = 11;  % rotated around vertical axis away from the nose, deg
    end
    
    properties            % variable depending on accommodation
        PupilD   = 3      % 3 mm pupil diameter by default
        LensD    = 10.6318 % value for 0D accommodation
                          % as measured by Hermans et al. (1990) and focus of 2 micrometers
        Lens1x   = -9.70  % depends on accommodation
        Lens1R   = 11.51  % front lens 0D radius (30-year old Goncharov 2009), 11.25 (Escudero-Sanz et al. 1999)
        Lens2x   = -5.70  % depends on accommodation
        Lens2R   = -6.01  % back lens 0D radius (30-year old Goncharov 2009), -5.70 (Escudero-Sanz et al. 1999)
        image = [];
    end
    
    methods
        
        function self = Lens( self, lens_diameter )
            % Create eye lens
            self.LensD = lens_diameter;            
            D = lens_diameter;
            
            % hyperbolic front surface
            V = self.LensV1;
            k = self.Lens1k;
            a = 1 + k;        
            A = 5 * D^6 * pi^2 - 1024 * a * V^2;
            d = a^3 * D^6 * ( 5 * D^12 * pi^4 + 2112 * a * D^6 * pi^2 * V^2 + 1572864 * a^2 * V^4 );
            B = ( -360 * a^2 * D^6 * pi^2 * V - 32768 * a^3 * V^3 + 5 * pi * sqrt( d ) )^(1/3);
            C = 10 * pi * D^2;
            % below are the 3 real roots of the cubic equation for the surface height, the 3d root gives the height
            %h1 = ( 32 * V + A / B - B / a  ) / C;  
            %h2 = ( 64 * V + ( 1 + 1i * sqrt(3) ) * (-A) / B + ( 1 - 1i * sqrt(3) ) * B / a ) / ( 2 * C );
            h3 = ( 64 * V + 1i * ( 1i + sqrt(3) ) * A / B + ( 1 + 1i * sqrt(3) ) * B / a ) / ( 2 * C );
            hf = real( h3 );  % remove zero imaginary part
            self.Lens1R = a * hf / 2 + D^2 / ( 8 * hf );
            
            % parabolic back surface
            V = self.LensV2;
            self.Lens2R = -pi * D^4 / ( 64 * V );
            hb = 8 * V / ( pi * D^2 );
            
            % update the lens eye element shape if the lens already exists
            if self.cnt > 3
                % find the current lens position (center of the equatorial circle)
                lv = self.elem{ 5 }.r - self.elem{ 4 }.r; % vector from the lens front to the lens back
                lc = self.elem{ 4 }.r + ( self.LensMx - self.Lens1x ) * lv / norm( lv );
                
                % create, orient, and position the new front surface
                lens = Lens( [ -hf 0 0 ], self.LensD, self.Lens1R, self.Lens1k, { 'aqueous' 'lens' } );
                lens.rotate( self.elem{ 4 }.rotax, self.elem{ 4 }.rotang ); % rotate the surface normal
                lens.r = rodrigues_rot( lens.r, self.elem{ 4 }.rotax, self.elem{ 4 }.rotang ); % rotate the surface position
                lens.r = lens.r + lc ;
                self.elem{ 4 } = lens;

                % create, orient, and position the new back surface
                lens = Lens( [ hb 0 0 ], self.LensD, self.Lens2R, self.Lens2k, { 'lens' 'vitreous' } );
                lens.rotate( self.elem{ 5 }.rotax, self.elem{ 5 }.rotang );
                lens.r = rodrigues_rot( lens.r, self.elem{ 4 }.rotax, self.elem{ 4 }.rotang );
                lens.r = lens.r + lc ;
                self.elem{ 5 } = lens;
            end
            self.Lens1x = self.LensMx - hf; % front surface apex x
            self.Lens2x = self.LensMx + hb; % back surface apex x
        end
        
        function self = Eye( lens_diameter, pupil_diameter )
            % Eye constructor function
            if nargin == 0 || ( nargin > 0 && isempty( lens_diameter ) )
                lens_diameter = self.LensD; % 0D accommodation by default
            end
            if nargin > 1
                self.PupilD = pupil_diameter;
            end
            
            % compound eye elements
            
            % cornea
            cornea1 = Lens( [ self.Cornea1x 0 0 ], self.Cornea1D, self.Cornea1R, self.Cornea1k, { 'air' 'cornea' } );
            self.cnt = self.cnt + 1;
            self.elem{ self.cnt } = cornea1;
            cornea2 = Lens( [ self.Cornea2x 0 0 ], self.Cornea2D, self.Cornea2R, self.Cornea2k, { 'cornea' 'aqueous' } );
            self.cnt = self.cnt + 1;
            self.elem{ self.cnt } = cornea2;
            
            % pupil
            % pupil = Aperture( [ self.Pupilx 0 0 ], [ self.PupilD, self.PupilDout, 1000 ] ); % the third diameter is the one actually used for tracing, the second just to draw
            pupil = Lens( [ self.Pupilx 0 0 ], [ self.PupilD, -1.77*self.RetinaR ], -self.RetinaR, self.Retinak, { 'cornea' 'soot' } ); % pupil extends around the retina to simulate the opaque choroid
            self.cnt = self.cnt + 1;
            self.elem{ self.cnt } = pupil;
            
            % lens
            % caluclate lens dimensions from its diameter assuming constant volume
            self.Lens( lens_diameter );
            
            lens1 = Lens( [ self.Lens1x 0 0 ], self.LensD, self.Lens1R, self.Lens1k, { 'aqueous' 'lens' } );
            self.cnt = self.cnt + 1;
            self.elem{ self.cnt } = lens1;
            lens2 = Lens( [ self.Lens2x 0 0 ], self.LensD, self.Lens2R, self.Lens2k, { 'lens' 'vitreous' } );
            self.cnt = self.cnt + 1;
            self.elem{ self.cnt } = lens2;
            
            % retina
            retina = Retina( [ self.Retinax 0 0 ], self.RetinaR, self.Retinak );
            % retina.rotate( [ 0 0 1 ], self.RetinaRotz * pi / 180 );
            self.cnt = self.cnt + 1;            
            self.elem{ self.cnt } = retina;
            self.image = retina.image;
        end
                
        function V = eye_lens_vol( self, D )
            % calculate eye lens volume given its surface parameters at 0D, this
            % is a utility function to find D such that total V = 160 mm^2
            x = self.Lens1R / ( 1 + self.Lens1k );
            h = x + sqrt( x^2 - D^2 / ( 4 * ( 1 + self.Lens1k ) ) );
            V1 = pi * D^2 / 8 * h * ( 1 - h / ( 6 * self.Lens1R / ( 1 + self.Lens1k ) + 3 * h ) );
            h = abs( D^2 /( 8 * self.Lens2R ) );
            V2 = pi / 8 * D^2 * h; % parabolic volume
            V = V1 + V2;
        end

    end
    
end

