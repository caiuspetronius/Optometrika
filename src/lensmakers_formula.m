function [ P, nrefr ] = lensmakers_formula( Rf, Rb, wavelength, glass  )
% LENSMAKERS_FORMULA calculates refractive index and approximate the lens power near 
% its axis using the Lensmaker's formula
%
% INPUT:
%    Rf - front lens surface curvature radius, meters
%    Rb - back lens surface curvature radius, meters
%    wavelenght - of lihght, meters
%    glass - optical medium of the lens, see refrindx.m for the list
%
% OUTPUT:
%    P - lens power in diopters
%    nrefr - lens refractive index for the given wavelength
%
% Copyright: Yury Petrov, 2016
%

nrefr = refrindx( wavelength, glass );
P = 1000 * ( nrefr - 1 ) * ( 1 / Rf - 1 / Rb + ( nrefr - 1 ) * h / ( nrefr * Rf * Rb ) ); % Diopters

end

