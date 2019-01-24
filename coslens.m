function x = coslens( y, z, args, flag )
% To trace through a general lens profile one needs to create a function,
% which takes two arguments ( y, z ) defining a position in the lens
% plane, an arbitrary number of additional arguments provided in the cell 
% array args, and, finally, a flag argument. On flag == 0 the function 
% should return the lens height x for the given position ( y, z ). Otherwise,
% the function should return the lens normal at this position. By convention, 
% the normal should point along the x-axis, i.e. in the same general 
% direction as the traced ray.
%
% COSLENS defines a radially symmetric cosine profile used in example4.
% the first two input arguments are coordinates in the lens plane, the
% third argument is a cell array holding the lens height argv{1}, and the
% cosine period argv{2}.

r = sqrt( y.^2 + z.^2 );
if flag == 0
    x = args{1} * ( 1 - cos( 2 * pi / args{2} * r ) );
else
    c = 1 ./ sqrt( 1 + ( 2 * pi * args{1} / args{2} .* sin( 2 * pi / args{2} * r ) ).^2 );
    s = sqrt( 1 - c.^2 );
    th = atan2( z, y );
    x = -sign( args{ 1 } ) * [ -sign( args{ 1 } ) * c, s .* cos( th ), s .* sin( th ) ];
end
end
