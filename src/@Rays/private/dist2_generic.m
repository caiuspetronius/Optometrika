function d = dist2_generic( l, r0, e, surf )  rend = r0 + l * e; % the ray's end     d = ( rend( :, 1 ) - surf.eval( rend( :, 2 ), rend( :, 3 ),sign(l*e(1))) ).^2;     %add information about direction in beam direction of xend