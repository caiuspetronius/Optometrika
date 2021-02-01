function f = exGaussian2D_fill( X, Y, Mu, Std, Lambda, angle )
% EXGAUSSIAN implements exponentially modified Gaussian distribution 
% http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
% extending it to 2D
% 
% Copyright: Yury Petrov, 2016
%

R = [ X(:) Y(:) ];
Rt = R * [ cos( -angle ), sin( -angle ); -sin( -angle ), cos( -angle ) ];
if Lambda( 2 ) < 0
    Rt( :, 2 ) = -Rt( :, 2 );
    Lambda( 2 ) = - Lambda( 2 );
end
f = exGaussian( Rt( :, 1 ), Mu( 1 ), Std( 1 ), Lambda( 1 ) ) .* exGaussian( Rt( :, 2 ), Mu( 2 ), Std( 2 ), Lambda( 2 ) );
f = f / sum( f );
f = reshape( f, size( X ) );