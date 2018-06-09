function [ data ] = UpdateData( data,volume )
%UPDATEDATA updates data structure variable to avoid repetitive and
%expensive computational evaluations. Computes object-function, exit-wave
%and far-field propagated wave from a 3D volume of -delta+ibeta
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018


%Complex line integrals over the sample volume
frames=Radon3D(volume,data);
%Compute complex object transmissivity function:
data.O=exp(1i.*data.k*frames);
%Multiply by probe function to get exit-wave:
data.phi0=repmat(data.probe,[1,1,size(frames,3)]).*(data.O);
%Propagate exit-wave in far-field:
data.PHI0=Ft2(data.phi0);

end

