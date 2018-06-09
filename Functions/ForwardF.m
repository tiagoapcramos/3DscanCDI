function [ Fx ] = ForwardF(volume,data )
%FORWARDF - function to perform forward model equivalent to a linear
%operator defined as F: N -> I =abs(fft2(probe*object))^2
% volume              - 3D volume / Tomographic Reconstruction - can be
%                       complex
% parameters          - Matrix with projection parameters for all
%                       projections, in the form [Theta, u, v, alpha, beta]
%                       (:,1) - projection angles (degrees)
%                       (:,2) - horizontal granslations (Pixels)
%                       (:,3) - vertical translations (Pixels)
%                       (:,4) - alpha (Degrees)
%                       (:,5) - beta (Degrees)
% data.sp             - Frame Size= probe array size [Pixels]
% data.max_memory     - Maximum available memory in GPU
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

if ~isfield(data,'k')
    try
        data.k = (2*pi*data.pixsize)/data.lambda;
    catch
        display('No k defined. using k=1')
        data.k = 1;
    end
end
if ~isfield(data,'probe')
    data.probe=false;
end

%Apply Radon transform to complex array 'volume'
frames=Radon3D(volume,data);
%Calculate Object Function
O=exp(1i.*data.k*frames);
%Multiply by probe function and propagate in far-field
if data.probe==false
    Fx=Ft2(O);
else
    Fx=Ft2(repmat(data.probe,[1,1,size(frames,3)]).*O);
end
%Compute Intensity of diffraction patterns:
Fx=abs(Fx).^2;


end

