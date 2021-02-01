function [rinter] = sphere_intersection( r_in, e, surf )
%
% returns intersection vector with full sphere surface

n1=e(:,1);
n2=e(:,2);
n3=e(:,3);
p1=r_in(:,1);
p2=r_in(:,2);
p3=r_in(:,3);
ps1=0;
ps2=0;
ps3=0;
rad_sph=surf.D(1)/2;

t(:,1)=  ((- n1.^2.*p2.^2 + 2.*n1.^2.*p2.*ps2 - n1.^2.*p3.^2 + 2.*n1.^2.*p3.*ps3 - n1.^2.*ps2.^2 - n1.^2.*ps3.^2 + n1.^2.*rad_sph.^2 + 2.*n1.*n2.*p1.*p2 - 2.*n1.*n2.*p1.*ps2 - 2.*n1.*n2.*p2.*ps1 + 2.*n1.*n2.*ps1.*ps2 + 2.*n1.*n3.*p1.*p3 - 2.*n1.*n3.*p1.*ps3 - 2.*n1.*n3.*p3.*ps1 + 2.*n1.*n3.*ps1.*ps3 - n2.^2.*p1.^2 + 2.*n2.^2.*p1.*ps1 - n2.^2.*p3.^2 + 2.*n2.^2.*p3.*ps3 - n2.^2.*ps1.^2 - n2.^2.*ps3.^2 + n2.^2.*rad_sph.^2 + 2.*n2.*n3.*p2.*p3 - 2.*n2.*n3.*p2.*ps3 - 2.*n2.*n3.*p3.*ps2 + 2.*n2.*n3.*ps2.*ps3 - n3.^2.*p1.^2 + 2.*n3.^2.*p1.*ps1 - n3.^2.*p2.^2 + 2.*n3.^2.*p2.*ps2 - n3.^2.*ps1.^2 - n3.^2.*ps2.^2 + n3.^2.*rad_sph.^2).^(1./2) - n1.*p1 - n2.*p2 - n3.*p3 + n1.*ps1 + n2.*ps2 + n3.*ps3)./(n1.^2 + n2.^2 + n3.^2);
t(:,2)= -((- n1.^2.*p2.^2 + 2.*n1.^2.*p2.*ps2 - n1.^2.*p3.^2 + 2.*n1.^2.*p3.*ps3 - n1.^2.*ps2.^2 - n1.^2.*ps3.^2 + n1.^2.*rad_sph.^2 + 2.*n1.*n2.*p1.*p2 - 2.*n1.*n2.*p1.*ps2 - 2.*n1.*n2.*p2.*ps1 + 2.*n1.*n2.*ps1.*ps2 + 2.*n1.*n3.*p1.*p3 - 2.*n1.*n3.*p1.*ps3 - 2.*n1.*n3.*p3.*ps1 + 2.*n1.*n3.*ps1.*ps3 - n2.^2.*p1.^2 + 2.*n2.^2.*p1.*ps1 - n2.^2.*p3.^2 + 2.*n2.^2.*p3.*ps3 - n2.^2.*ps1.^2 - n2.^2.*ps3.^2 + n2.^2.*rad_sph.^2 + 2.*n2.*n3.*p2.*p3 - 2.*n2.*n3.*p2.*ps3 - 2.*n2.*n3.*p3.*ps2 + 2.*n2.*n3.*ps2.*ps3 - n3.^2.*p1.^2 + 2.*n3.^2.*p1.*ps1 - n3.^2.*p2.^2 + 2.*n3.^2.*p2.*ps2 - n3.^2.*ps1.^2 - n3.^2.*ps2.^2 + n3.^2.*rad_sph.^2).^(1./2) + n1.*p1 + n2.*p2 + n3.*p3 - n1.*ps1 - n2.*ps2 - n3.*ps3)./(n1.^2 + n2.^2 + n3.^2);

% Do not return 0 value
% Return minimum positive distance

find( t>eps);
t(t<eps)=NaN;
d=min(t,[],2);

% form the intersection vector
rinter = r_in + repmat( d, 1, 3 ) .* e;
                        

end