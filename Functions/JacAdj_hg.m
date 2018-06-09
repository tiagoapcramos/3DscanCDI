function [ Ja_hg ] = JacAdj_hg(volume,data,hg )
%JACADJ_HG - Applies the Jacobiant Adjoint operator to the increment hg
% volume              - 3D volume / Tomographic Reconstruction - can be
%                       complex
% parameters          - Matrix with projection parameters for all
%                       projections, in the form [Theta, u, v, alpha, beta]
%                       (:,1) - projection angles (degrees)
%                       (:,2) - horizontal granslations (Pixels)
%                       (:,3) - vertical translations (Pixels)
%                       (:,4) - alpha (Degrees)
%                       (:,5) - beta (Degrees)
% data.sp             - Frame Size= probe array size
% data.max_memory     - Maximum available memory in GPU
% Ja_hg               - Result after applying the linear operator J' to hg
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
    data.probe=ones(data.sp);
end

if ~isfield(data,'phi0')
    % First apply Radon transform to complex array 'volume'
    frames=Radon3D(volume,data);
    %Compute complex object transmissivity function
    O=exp(1i.*data.k*frames);
    %Multiply by probe function
    phi0=repmat(data.probe,[1,1,size(frames,3)]).*(O);
else
    phi0=data.phi0;
end
    %Free-space Propagation
if ~isfield(data,'PHI0')
    PHI0=Ft2(phi0);
else
    PHI0=data.PHI0;
end
F1=conj(1i.*data.k.*phi0);
if isfield(data,'noise_model') 
    if strcmp(data.noise_model,'poisson')
      hg=hg./(sqrt(data.Imeas+1));
    end
end
F2=Ift2(PHI0.*real(hg));
F3=F1.*F2;
Ja_hg=2*BackRadon3D(F3,data);
Ja_hg(isnan(Ja_hg))=0;

end

