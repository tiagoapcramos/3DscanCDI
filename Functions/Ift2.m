function [ X ] = Ift2( F,correction)
%IFT2 - Performs a 2D Discrete Inverse Fourier Transform on the stack of 
%arrays F. Takes into account shifting of the zero-frequency component to
%the center of the spectrum, and propper scalling according to Parseval's 
%theorem.
% input:
% F          -   Matrix to apply inverse transformation. If size(F,3)>1 F 
%                is interpreted as a stack of matrices in its 3rd direction  
% correction -  (optional): if 'symmetric' causes Ift2 to treat F as
%               conjugate symmetric. if 'real' or 'abs' take the real part
%               or absolute value of the resulting transformation. default:
%               'none'
% output:
% X     -   2D Discrete Inverse Fourier Transform of F.
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

[Ny,Nx,~]=size(F);

if nargin<2 ||isempty(correction),correction='none';end

switch correction
    case 'symmetric'
    X=fftshift(fftshift(ifft2(ifftshift(ifftshift(F,1),2),'symmetric'),1),2);
    case 'real'
    X=real(fftshift(fftshift(ifft2(ifftshift(ifftshift(F,1),2)),1),2));
    case 'abs'
    X=abs(fftshift(fftshift(ifft2(ifftshift(ifftshift(F,1),2)),1),2));
    case 'none'
    X=fftshift(fftshift(ifft2(ifftshift(ifftshift(F,1),2)),1),2);
end

X=X*sqrt(Nx*Ny);


end

