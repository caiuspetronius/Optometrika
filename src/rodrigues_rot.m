function v_rot = rodrigues_rot( v, k, theta )
% RODRIQUES_ROT - Rotates array of 3D vectors by an angle theta about vector k.
% Direction is determined by the right-hand (screw) rule.
%
% Syntax:  v_rot = rodrigues_rot( v, k, theta )
%
% INPUT:
%    v     - n x 3 array of vectors to be rotated
%    k     - 1 x 3 rotation axis
%    theta - rotation angle, radians
% OUTPUT:
%    v_rot - Array of rotated vectors.
%
% Copyright: Yury Petrov, 2016
%

[ m, n ] = size( v );
if n ~= 3
    error( 'Input vector is not 3 dimensional!' );
end
if size( k, 2 ) ~= 3
    error( 'Rotation axis is not 3 dimensional!' );
end

k = k / norm( k ); % normalize rotation axis
k = repmat( k, m, 1 );
v_rot = v .* cos( theta ) + cross( k, v, 2 ) .* sin( theta ) + k .* repmat( dot( k, v, 2 ), 1, 3 ) .* ( 1 - cos( theta ) );
end