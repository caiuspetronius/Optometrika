function example18()
% Analyze complete geometrical OPL error and offset error 
% for a beam that passes through two non-wedged windows
% that are tilted +45° and -45 degree.
% A. Schultze 2020-10-05


%Add src to PATH
startup;

% create a container for optical elements (Bench class)
bench = Bench;
ref_bench = Bench;
% add optical elements in the order they are encountered by light rays

medium = 'vacuum';



bs1a = Plane( [ 10 0 0 ] , 100*sqrt(2), { medium  'bk7' } );
bs1b = Plane( [ 10+20 0 0 ] , 100*sqrt(2), {'bk7' medium} );
bs2a = Plane( [ 60 0 0 ] , 100*sqrt(2), { medium  'bk7' } );
bs2b = Plane( [ 60+20 0 0 ] , 100*sqrt(2), {'bk7' medium});

icase = '1bs untilted';
switch(icase)
    case '1bs untilted'
        bench.append( bs1a);
        bench.append( bs1b);
    case '1bs tilted'    
        bs1a.rotate([0 0 1], 1*pi/4);
        bs1b.rotate([0 0 1], 1*pi/4); 
        bench.append( bs1a);
        bench.append( bs1b);
    case '2bs tilted' 
        bs1a.rotate([0 0 1], 1*pi/4);
        bs1b.rotate([0 0 1], 1*pi/4);
        bs2a.rotate([0 0 1], 1*-pi/4);
        bs2b.rotate([0 0 1], 1*-pi/4);
        bench.append( bs1a);
        bench.append( bs1b);
        bench.append( bs2a);
        bench.append( bs2b);
end


screenwf = ScreenGeneric([ 100 0 0 ], 100, 100, 512, 512,'wf' );
bench.append( screenwf );

ref_bench.append( screenwf );

%variation in position
rays_in = Rays( 5, 'source', [ 1 0 0 ], [ 1 0 0 ], 10e-6, 'square',medium,1064e-9,[],0);
rays_out = bench.trace_recursive( rays_in );
rays2_out = ref_bench.trace( rays_in );
bench.draw(rays_out, 'lines');

%variation in angle




figure();
opl = rays_out(end).opl;
opl_ref = rays2_out(end).opl;

phi_y = asin(rays_in.n(:,2)./rays_in.n(:,1));
phi_z = asin(rays_in.n(:,3)./rays_in.n(:,1));
subplot(2,2,1);
plot3(phi_y,phi_z,(opl-opl_ref)-mean(opl-opl_ref),'.');
xlabel('Beam Angle Y (rad)');
ylabel('Beam Angle Z (rad)');
zlabel('Optical Path Length (mm)');
title('Delta OPL until screen (referenced)');

subplot(2,2,2);
plot3(phi_y,phi_z,opl_ref-mean(opl_ref),'.');
xlabel('Beam Angle Y (rad)');
ylabel('Beam Angle Z (rad)');
zlabel('Optical Path Length (mm)');
title('Delta OPL until screen (reference w/o glass)');



off = vecnorm(rays_out(end).r(:,2:3)');
off_ref = vecnorm(rays2_out(end).r(:,2:3)');
phi_y = asin(rays_in.n(:,2)./rays_in.n(:,1));
phi_z = asin(rays_in.n(:,3)./rays_in.n(:,1));
subplot(2,2,3);

plot3(phi_y,phi_z,off-off_ref,'.');
xlabel('Beam Angle Y (rad)');
ylabel('Beam Angle Z (rad)');
zlabel('Offset (mm)');
title('Beam Offset on screen ( referenced)');
subplot(2,2,4);
plot3(phi_y,phi_z,off_ref,'.');
%plot(theta,opl-mean(opl));
xlabel('Beam Angle Y (rad)');
ylabel('Beam Angle Z (rad)');
zlabel('Offset (mm)');
title('Beam Offset on screen (reference w/o glass)');


suptitle(['Example 18: ' icase]);