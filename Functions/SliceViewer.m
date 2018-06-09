function [ ] = SliceViewer(X,alph)
% SLICEVIEWER Makes image of sagital coronal and axial central slices of 
% tomographic volume X
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018


if nargin<2||isempty(alph),alph=0.5;end
rc=isreal(X(:));%Checks if volume is real or complex.
mrX=min(real(X(:)));
miX=max(imag(X(:)));
% X(real(X(:))==0)=NaN;
axi=squeeze(X(:,:,round(end/2)));
sagi=squeeze(X(:,round(end/2),:));
coron=squeeze(X(round(end/2),:,:));


[ny,nx,~]=size(X);
[x1,y1]=meshgrid(linspace(-nx/2,nx/2,nx),linspace(-ny/2,ny/2,ny));
if rc==0
   ha=tight_subplot(1,2,[0 0.1]);
   hold on
   axes(ha(1)),surface(x1,y1,x1*0,real(axi)/mrX,'FaceColor',...
       'texturemap','EdgeColor','none','CDataMapping','scaled');
   axes(ha(1)),surface(x1,x1*0,y1,real(sagi)'/mrX,'FaceColor',...
       'texturemap','EdgeColor','none','CDataMapping','scaled');
   axes(ha(1)),surface(x1*0,x1,y1,real(coron)'/mrX,'FaceColor',...
       'texturemap','EdgeColor','none','CDataMapping','scaled');
   title('Real Part')
   axes(ha(2)),surface(x1,y1,x1*0,imag(axi)/miX,'FaceColor',...
       'texturemap','EdgeColor','none','CDataMapping','scaled');
   axes(ha(2)),surface(x1,x1*0,y1,imag(sagi)'/miX,'FaceColor',...
       'texturemap','EdgeColor','none','CDataMapping','scaled');
   axes(ha(2)),surface(x1*0,x1,y1,imag(coron)'/miX,'FaceColor',...
       'texturemap','EdgeColor','none','CDataMapping','scaled');
   title('Imaginary part')

   makeitpretty
   for k=1:2
       axes(ha(k)),
       makeitpretty
       view([32 31])
       axis equal tight,
       alpha(alph)
       set(gca,'Xtick',[])
       set(gca,'Ytick',[])
       set(gca,'Ztick',[])
       ax=gca;
       ax.BoxStyle='full';
   end
else
   hold on
   surface(x1,y1,x1*0,real(axi)/mrX,'FaceColor','texturemap',...
       'EdgeColor','none','CDataMapping','scaled');
   surface(x1,x1*0,y1,real(sagi)'/mrX,'FaceColor','texturemap',...
       'EdgeColor','none','CDataMapping','scaled');
   surface(x1*0,x1,y1,real(coron)'/mrX,'FaceColor','texturemap',...
       'EdgeColor','none','CDataMapping','scaled');
   title('Real Part')
   makeitpretty
   view([32 31])
   axis equal tight,
   alpha(alph)
   set(gca,'Xtick',[])
   set(gca,'Ytick',[])
   set(gca,'Ztick',[])
   ax=gca;
   ax.BoxStyle='full';
end

end

