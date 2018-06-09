function [probe] = GaussianProbe(sp,c1,c2,c3)
%GAUSSIANPROBE generates a complex array of a probe function to be used for
%3D ptychogrpahy data generation
% INPUT:
%   sp -    Size of probe function [rows,collumns];
%   c1 -    shape constant - variance/FWHM
%   c2 -    phase constant - Phase Variation/Gradient coefficient
%   c3 -    intensity constant - Maximum probe amplitude
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

if nargin<2||isempty(c1),c1=5;end
if nargin<3||isempty(c2),c2=100;end
if nargin<4||isempty(c3),c3=1e-3;end

if numel(sp)==1
    sp=[sp,sp];
end
probe=fspecial('gaussian',max(sp),max(sp)/c1)+eps;
probe=probe(1:sp(1),1:sp(2));
probe=sqrt(probe);
probe=probe/max(probe(:));
probe=probe.*exp(1i.*c2*(probe))*c3;


end

