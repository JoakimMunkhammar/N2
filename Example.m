%        DownscalingTest.m
%-----------------------------------
%
% This code runs the N2Downscaling
% code with an example input on hour
% resolution global horizontal irradiance
% clear-sky index and example
% output in minute resolution global,
% beam and diffuse horizontal irradiance
% clear-sky indices.

% An example mean value GHI clear-sky index
MeanGHICSI = [0.5 0.1 0.3 0.9 0.9 0.8 0.3];

% Using the N2Downsccaling model
[Global,Beam,Diffuse] = N2Downscaling(MeanGHICSI);