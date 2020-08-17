function focal = example17_wavefront()
% Measure the pathlength variation (OPL) of three different lens types.
% A. Schultze 2020-08-17


for i=1:3 %each lens type

% create a container for optical elements (Bench class)
bench = Bench;

% add optical elements in the order they are encountered by light rays

lenses= {'Plano Convex','Best Form','Aspheric'};

if i==1
    % 1) Plano Convex Standard Lens, 100mm #LA1509-C
    wd = 99.7;
    lens1 = Plane( [ wd 0 0 ]   , 25.4, { 'air' 'bk7' } ); % flat
    lens2 = Lens( [ wd+3.6 0 0 ], 25.4, -51.5, 0, { 'bk7' 'air' } ); % convex

elseif i==2
    % 2) Best Form Lens #Thorlabs LBF254-100, 100mm
    lens1 = Lens( [ 100 0 0 ]    , 25.4, 60.02, 0, { 'air' 'bk7' } ); 
    lens2 = Lens( [ 100+4.0 0 0 ], 25.4, -353.30, 0, { 'bk7' 'air' } ); 
elseif i==3
    % 3) Aspheric Lens #Thorlabs ASL10142-C, EFL=79.0mm
    wd = 74.84+2.1;
    lens1 = Plane( [ wd 0 0 ]   , 25.4, { 'air' 'fusedsilica' } ); % parabolic surface 
    lens2 = AsphericLens( [ wd+3.5 0 0 ], 25.4, -35.8256, -0.6291,-[0 1.439e-7], { 'fusedsilica' 'air' } ); 
end
 

bench.append( { lens1, lens2 } );

% screen
screen = Screen( [ 110 0 0 ], 25, 25, 512, 512 );
bench.append( screen );

screenwf = ScreenWavefront( [ 115 0 0 ], 20, 20, 100, 100,1064e-9 );
bench.append( screenwf );

%bench.rotate( [ 0 0 1 ], 0.15 );

% create some rays
nrays = 1000;

rays_in = Rays( nrays, 'source', [ 0 0 0 ], [ 1 0.0 0 ], 0.1, 'square','air',1064e-9);


tic;

rays_through = bench.trace( rays_in );    % repeat to get the min spread rays

% draw bench elements and draw rays as arrows

if i==1
    bench.draw( rays_through, 'lines' );  % display everything, the other draw option is 'lines'
end

% Wavefront plots
screenwf.plot3();
title(['Lens ' num2str(i) ': ' lenses{i}]);
%screenwf.surf();


end %for i=1:3 %each lens type

if(0)
%Show the "image" of the wavefront image
%This is a less elegant version of what is done above using hist2mean.
figure( 'Name', 'Wavefront', 'NumberTitle', 'Off' );
him=imshow( screenwf.image, [] );
set(him, 'AlphaData', 1-0.95*(screenwf.image==0))
title('Wavefront Pathlength in mm');
colorbar;
colormap default;
end

end

