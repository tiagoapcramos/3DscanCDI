function [ F ] = Ft2(X)
%FT2 - Performs a 2D Discrete Fourier Transform on the stack of arrays X.
%Takes into account shifting of the zero-frequency component to the center
%of the spectrum, and propper scalling according to Parseval's theorem
% input:
% X     -   Matrix to apply transformation. If size(X,3)>1 X is interpreted
%           as a stack of matrices in its 3rd direction  
% output:
% F     -   2D Discrete Fourier Transform of X.
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

[Ny,Nx,~]=size(X);
F=fftshift(fftshift(fft2(ifftshift(ifftshift(X,1),2)),1),2);
F=F/sqrt(Nx*Ny);

end

