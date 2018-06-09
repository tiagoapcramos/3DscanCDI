function [ mask ] = Mask_Circle(D,s )
%MASK_CIRCLE - Creates a 2D binarry array of size [s,s] Pixels with a
%circular mask of diameter D (Pixels).
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

if nargin<2||isempty(s)
    s=[D D];
end
if size(s,2)==1
    s=[s s];
end
[X,Y]=meshgrid(1:s(1),1:s(2));
cx=ceil(s(1)/2);
cy=ceil(s(2)/2);
circle=(X-cx).^2+(Y-cy).^2;
mask=circle<=(D/2)^2;

end

