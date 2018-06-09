function [f] = imloop( X )
%IMLOOP images real or complex arrays. If X is a real number a single image
%frame is generated. if X is complex thn 4 different images (real,
%imaginary, amplitude and phase) are displayed instead. If size(X,3)>1 then
%IMLOOP interprets X as a stack of images in its 3rd dimension.
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

ir=isreal(X); %is real
f=figure('color','w');
if ir
    for k=1:size(X,3)
        imagesc(X(:,:,k))
        axis equal tight
        title(num2str(k))
        pause(10/size(X,3))
    end
else
    for k=1:size(X,3)
        subplot(221),imagesc(real(X(:,:,k))),axis equal tight, title(['real',num2str(k)])
        subplot(222),imagesc(imag(X(:,:,k))),axis equal tight, title('imag')
        subplot(223),imagesc(abs(X(:,:,k))),axis equal tight, title('abs')
        subplot(224),imagesc(angle(X(:,:,k))),axis equal tight, title('angle')
        pause(10/size(X,3))
    end
end
    



end

