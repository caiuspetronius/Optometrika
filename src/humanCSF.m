function csf = humanCSF( f )
%
% This function approximates human contrast sensitivity function (CSF), 
% which is definced as the reciprocal of the threshold contrast. The
% analytical fit is from Watson AB and Ahumada AJ,Jr. "A standard model for 
% foveal detection of spatial contrast." J Vis. 2005; 5 (9): 717?740.
% The analytical form is given by the difference of hyperbolic secants 
% describing high and low frequency parts of CSF and defined by 5 parameters:
%
% CSF( f; f0,f1,a,p ) = gain * ( sech((f/f0)^p) ? a * sech(f/f1) )
%

gain = 373.08;
f0 = 4.1726;
f1 = 1.3625;
a = 0.8493;
p = 0.7786;
% beta = 2.4081; % summation over channels metric, not used in CSF fit
% sigma = 0.6273; % summation over periphery of the fovea, not used	

csf = gain * ( sech( ( f / f0 ).^p ) - a * sech( f / f1 ) );