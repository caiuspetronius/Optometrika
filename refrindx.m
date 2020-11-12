function [ n_refr, dens, hard ] = refrindx( wavelength, glass )
% REFRINDX refractive index calculations.
%
% INPUT:
%   wavelength - light wavelength in meters
%   glass - material string
%
% OUTPUT:
%   n_refr - refractive index for the wavelengths: 587.6, 486.1, 656.3 nm
%   dens - material density in g/cm^3
%   hard - abrasion resistance from 1 to 10 (Mohs hardness)
%
% Data from JML Optical Industries, Inc available at:
%   http://www.netacc.net/~jmlopt/transmission2.html
%   Or at:
%   http://www.jmlopt.com/level2/TechInfo/materialstable.aspx
%
%   The human eye data are from in Escudero-Sanz & Navarro, 
%   "Off-axis aberrations of a wide-angle schematic eye model", 
%   JOSA A, 16(8), 1881-1891 (1999).

persistent glass_names_sellmeier refr_consts_sellmeier

if isempty( wavelength )
    wavelength = 5876e-10; % green
end

n_refr = [];
dens = [];
hard = [];

if isempty( glass_names_sellmeier )
    qwe = '';
    fp = fopen( 'Sellmeier.glass.refr', 'r' );
    while ~feof( fp )
        qwe = str2mat( qwe, fgetl( fp ) );
    end
    glass_names_sellmeier = qwe( 2:end, 1:7 );
    refr_consts_sellmeier = str2num( qwe( 2:end, 8:end ) );
end

lambda_ref = [ 5876 4861 6563 ] * 1e-10;
lambda_eye = [4580 5430 5893 6328] * 1e-10;

if isempty( wavelength )
    wavelength = lambda_ref( 1 ); % green wavelength by default
end

ind = strmatch( glass, glass_names_sellmeier );
if ~isempty( ind )
    lambda = wavelength * 1e6; % change units to micrometer
    A = refr_consts_sellmeier( ind, 1:3 );
    B = refr_consts_sellmeier( ind, 4:6 );
    n_refr = sqrt( 1 + A(1) * lambda.^2 ./ ( lambda.^2 - B(1) ) + ...
        A(2) * lambda.^2 ./ ( lambda.^2 - B(2) ) + ...
        A(3) * lambda.^2 ./ ( lambda.^2 - B(3) ) );
else
    switch lower( glass )
        case {'vacuum' 'mirror' 'soot'}
            nref = [1];
            dens = 0;
            hard = 0;
        case 'air'
            nref = [1.00027717 1.00027951 1.00027625];
            dens = 0.001225;
            hard = 0;
        case 'water'
            nref = [ 1.3325 1.3356 1.3310 ];
            dens = 1.;
            hard = 0;
        case 'ice'
            nref = [ 1.3098 1.3137 1.3079 ];
            dens = 0.9167;
            hard = 1.5;
        case { 'quartz' 'SiO2' }
            nref = [ 1.4585 1.4631 1.4564 ];
            dens = 2.648;
            hard = 7;
        case 'diamond'
            nref = [ 2.4168 2.4305 2.4088 ];
            dens = 3.51;
            hard = 10;
        % plastics
        case { 'pmma', 'acrylic' }
            nref = [ 1.491 1.496 1.488 ];
            dens = 1.185;
            hard = 10;
        case { 'pc', 'polycarbonate' }
            nref = [ 1.5849 1.5994 1.5782 ];
            dens = 1.21;
            hard = 2;
        case { 'ps', 'polystyrene' }
            nref = [ 1.5917 1.6056 1.5853 ];
            dens = 1.06;
            hard = 4;
        case { 'nas-21' 'nas' }
            nref = [ 1.5714 1.5835 1.5669 ];
            dens = 1.09;
            hard = 6;
        case { 'optorez-1330' 'optorez' }
            nref = [ 1.5094 1.5163 1.5067 ];
            dens = 1.19;
        case { 'zeonex-e48r' 'zeonex' }
            nref = [ 1.5309 1.5376 1.5273 ];
            dens = 1.01;
            hard = 10;
        case 'optimas6000'
            nref = [1.5013 1.5061 1.4969];
            dens = 1.08;
            hard = 10;
        case { 'r-5000' 'polysulfone' 'psu' }
            nref = [ 1.675 1.696 1.660];
            dens = 1.24;
        % glasses
        case 'b270'
            nref = [1.5230 1.5292 1.5202];
            dens = 2.55;
        case 'bak1'
            nref = [1.5725 1.5794 1.5694];
            dens = 3.19;
        case 'bak2'
            nref = [1.5399 1.5462 1.5372];
        case 'bak4'
            nref = [1.56883 1.5759 1.56576];
            dens = 3.10;
        case 'balkn3'
            nref = [1.51849 1.52447 1.51586];
        case 'bk7'
            nref = [1.5168 1.5224 1.5143];
            dens = 2.51;
        case 'f2'
            nref = [1.6200 1.6320 1.6150];
            dens = 3.60;
        case 'f3'
            nref = [1.61293 1.62461 1.60806];
        case 'f4'
            nref = [1.6165 1.6284 1.6116];
        case 'fusedsilica'
            nref = [1.458 1.463 1.456];
            dens = 2.20;
        case 'k5'
            nref = [1.5224 1.5285 1.5198];
            dens = 2.59;
        case 'k7'
            nref = [1.51112 1.517 1.50854];
        case 'lasfn9'
            nref = [1.850 1.8689 1.8425];
            dens = 4.44;
        case 'lah71'
            nref = [1.8502 1.8689 1.8425];
        case 'n-bk7'
            nref = [1.5197 1.5239 1.5156];
        case 'pyrex'
            nref = [1.473 1.478 1.471];
        case 'sapphire'
            nref = [1.7682 1.7756 1.7649];
            dens = 3.97;
        case 'sf1'
            nref = [1.71736 1.73462 1.71031];
            dens = 3.66;
        case 'sf2'
            nref = [1.6476 1.6612 1.6421];
            dens = 3.86;
        case 'sf5'
            nref = [1.6727 1.6875 1.66661];
            dens = 4.07;
        case 'sf8'
            nref = [1.6889 1.7046 1.6825];
            dens = 4.22;
        case 'sf10'
            nref = [1.72825 1.74648 1.72085];
            dens = 4.28;
        case 'sf11'
            nref = [1.7847 1.8064 1.7759];
            dens = 5.41;
        case 'sf12'
            nref = [1.64831 1.66187 1.64271];
        case 'sf15'
            nref = [1.69895 1.71546 1.69221];
            dens = 2.92;
        case 'sf18'
            nref = [1.7215 1.7390 1.7143];
            dens = 4.49;
        case 'sf19'
            nref = [1.6668 1.6811 1.6609];
        case 'sf56'
            nref = [1.7847 1.8061 1.7760];
            dens = 3.28;
        case 'sk3'
            nref = [1.6088 1.6160 1.6056];
        case 'sk5'
            nref = [1.5891 1.5958 1.5861];
        case 'sk11'
            nref = [1.5638 1.5702 1.5610];
            dens = 3.08;
        case 'sk16'
            nref = [1.6204 1.6275 1.6172];
        case 'ssk2'
            nref = [1.6223 1.63048 1.61878];
        case 'ssk4a'
            nref = [1.61765 1.62547 1.61427];
        case 'ssk51'
            nref = [1.60361 1.61147 1.60022];
        case 'zk5'
            nref = [1.53375 1.54049 1.53084];
        % eye
        case 'cornea'
            nref = [1.3828 1.3777 1.376 1.3747];
        case 'aqueous'
            nref = [1.3445 1.3391 1.3374 1.336];
        case 'lens'
            nref = [1.4292 1.4222 1.42 1.4183];
        case 'vitreous'
            nref = [1.3428 1.3377 1.336 1.3347];
        otherwise
            error( [ 'No values for glass absorption for: ', glass ] );
    end
    
    if size( nref ) == 1
        n_refr=nref*ones(size(wavelength));
    elseif strcmp( glass, 'cornea' ) || ...
            strcmp( glass, 'aqueous' ) || ...
            strcmp( glass, 'lens' ) || ...
            strcmp( glass, 'vitreous' )
        [ abc, ~, Mu ] = polyfit( 1./lambda_eye.^2, nref, 2 );
         n_refr = polyval( abc, 1./wavelength.^2, [], Mu );
    else
        [ abc, ~, Mu ] = polyfit( 1./lambda_ref.^2, nref, 2 );
         n_refr = polyval( abc, 1./wavelength.^2, [], Mu );
    end
   
end

end
