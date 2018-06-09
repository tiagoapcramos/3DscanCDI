function [ a ] = Save_Gif(k,file_name,loops )
if nargin<3||isempty(loops)
    loops=inf;
end
%drawnow

f=getframe(gcf);
im=frame2im(f);
[imind,cm]=rgb2ind(im,256);
if k==1
    imwrite(imind,cm,file_name,'gif', 'Loopcount',loops);
else
    imwrite(imind,cm,file_name,'gif','WriteMode','append');
end
a=k;
end

