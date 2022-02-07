function meas=example19()
% Demonstrate accuracy of the numerical solver (used for generic and aspheric lenses) 
% and show performance comparison for a flat surface for default
% (quadratic) mode and higher precision mode.
% Limitation of both is eps of the solver.


% REMARK: For reasonable results you should change the initial condition
% of the numeric solver to a random value
% Rays.m:501 -#  [ d, fval ] = fminunc( a_fun, 20, ...
% to         +#  [ d, fval ] = fminunc( a_fun, rand, ...

% A. Schultze 2021-01-28

%Add src to PATH
startup;

bench = Bench;

% front lens surface FLAT
lens1 = AsphericLens( [ 0 0 0 ], 31, 1e16, 0, [ 0 0 0 ], { 'air' 'N-BAK 4' } );
bench.append( lens1 );
tic;
nrays = 500;


%% Default Mode (Fast)
rays_in = Rays( nrays, 'collimated', [ -10 0 0 ], [ 1 0 0 ], 25, 'hexagonal','air' );
rays_in.setget_solver_flag('default');

rays_through = bench.trace( rays_in );    % repeat to get the min spread rays
toc;
% draw bench elements and draw rays as arrows
% bench.draw( rays_through );  % display everything, the other draw option is 'lines'

figure();
plot3(rays_through(2).r(:,2),rays_through(2).r(:,3),rays_through(2).r(:,1),'.')
xlabel('x (mm)');ylabel('y (mm)');zlabel('z-gpl (mm)');
legend(' default','location','SouthEast');view([90 0]);


disp('Error level for default quadratic minimisation (m)');
disp((max(rays_through(2).r(:,1))-min(rays_through(2).r(:,1)))/1e3);

%% High Precision Mode
rays_in.setget_solver_flag('precise');

rays_through = bench.trace( rays_in );    % repeat to get the min spread rays
toc;
figure();
plot3(rays_through(2).r(:,2),rays_through(2).r(:,3),rays_through(2).r(:,1),'.')
xlabel('x (mm)');ylabel('y (mm)');zlabel('z-gpl (mm)');
legend('precise','location','SouthEast');view([90 0]);

disp('Error level for precise quadratic minimisation (m)');
disp((max(rays_through(2).r(:,1))-min(rays_through(2).r(:,1)))/1e3);

%% High Precision Mode, Absolute Target Function
rays_in.setget_solver_flag('precise_abs');

rays_through = bench.trace( rays_in );    % repeat to get the min spread rays
toc;
figure();
plot3(rays_through(2).r(:,2),rays_through(2).r(:,3),rays_through(2).r(:,1),'.')
xlabel('x (mm)');ylabel('y (mm)');zlabel('z-gpl (mm)');
legend('precise abs','location','SouthEast');view([90 0]);

disp('Error level for precise absolute minimisation (m)');
disp((max(rays_through(2).r(:,1))-min(rays_through(2).r(:,1)))/1e3);