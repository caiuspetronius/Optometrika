%example21

%generate arbitrary function surface profiles
startup;

close all;


tic;

for i=1:3
    bench = Bench;
    switch i
        case 1
            % a flat mirror with a step profile
            mirror = SurfaceGeneric([10 00 0],[20 20],@(x,y) 0*x+0*y,{ 'air' 'mirror' } );
            prof = 2/3*[1 1 1; 1 1 1; 0 0 0; 0 0 0]; %step profile (45 deg)
            mirror.profile_set(prof);
            mirror.plot();
            suptitle('Step Surface 1');
        
        case 2
            % a flat mirror with a step profile
            mirror = SurfaceGeneric([10 00 0],[20 20],@(x,y) 0*x+0*y,{ 'air' 'mirror' } );
            prof = 0.2/3*[1 1 1; 1 1 1; 0 0 0; 0 0 0]; %step profile (0.45 deg)
            mirror.profile_set(prof);
            mirror.plot();
            suptitle('Step Surface 2');
        case 3
            %a complex mirror with a letter "E"
            mirror = SurfaceGeneric([10 0 0],[20 20], @(X,Y) 0.1*X .* exp(-X.^2 - Y.^2),{ 'air' 'mirror' } );
            prof = 0.03*[0 0 0 0;0 0 0 0; 0 1 1 0; 0 1 1 0; 0 1 0 0; 0 1 0 0; 0 1 1 0; 0 1 1 0;0 1 0 0;0 1 0 0; 0 1 1 0; 0 1 1 0; 0 0 0 0;0 0 0 0]; %letter E
            mirror.profile_set(prof);
            mirror.plot();
            suptitle('Surface 2 ("E")');
            
        case 4
            % sphere in cartesian coordinates
            r=20;
            fun_sph = @(x,y) r*sqrt(1-x.^2-y.^2);
            mirror = SurfaceGeneric([30 00 0],[40 40],fun_sph,{ 'air' 'mirror' } );
            mirror.plot();
            suptitle('Sphere');
            
    end %switch
    
    
    bench.append( mirror );
    screenwf = ScreenGeneric( [ 0 0 0 ], 20, 20, 100, 100,'wf' ); %other options are 'opl' 'wf' 'tilt'
    screenwf.rotate([0 1 0],pi);
    bench.append( screenwf );
    
    % create collimated rays
    nrays = 30;
    rays_in = Rays( nrays, 'collimated', [ 0 0 0 ], [ 1 0 0 ], 16, 'hexagonal','air' );
    toc;
    rays_through = bench.trace( rays_in );    % repeat to get the min spread rays
    toc;
    bench.draw( rays_through, 'lines' );  % display everything, scale arrow length 2x
    
    figure();
    screenwf.plot3();
end %for
