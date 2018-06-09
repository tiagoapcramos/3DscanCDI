function [ phantomt ] = newphantom2( N,M,x )
%NEWPHANTOM2 - generates a 3D volume phantom to be used in 3D ptychography
%simulations. The generated volume resembles a slightly conical capillary
%tube containing a sphere in its center (with an internal gradient) and a 
%number of cubes randomly oriented.
%Input:
% N     -   size of the simulated phantom [Pixels] (1st and 2nd dimension / width)
% M     -   size of the simulated phantom [Pixels] (3rd dimension / height)
% x     -   number of cubes inside the sample
% -----------------
%Output:
% phantomt -    3D simulated phantom of size [N,N,M]
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018
 
if nargin<2 || isempty(M)
    M=N;
end

phantomt=zeros(N,N,M);
[X,Y,Z]=meshgrid(1:N,1:N,1:M);
xc=N/2;
yc=xc;
zc=M/2;

%Design the outter capillary
r_out=linspace(0.42*N,0.38*N,M);
thickness=ceil(N/40);
r_in=r_out-thickness;
array2D=(X-xc).^2+(Y-yc).^2;
for k=1:M
    phantomt(:,:,k)=0.1*double(array2D(:,:,k)>=r_in(k).^2 & array2D(:,:,k)<=r_out(k).^2);
end

%Design the inner sphere
sphere_radius=5*thickness;
array3D=(X-xc).^2+(Y-yc).^2+(Z-zc).^2;

for k=floor(zc-sphere_radius):ceil(zc+sphere_radius)
    phantomt(:,:,k)=phantomt(:,:,k)+1/(sphere_radius^2)*array3D(:,:,k).*double(array3D(:,:,k)<=sphere_radius^2);
end


%Design first cube
cube_side=sphere_radius/2;
xmin=xc-cube_side/2;
ymin=yc-cube_side/2;
zmin=zc-cube_side/2;
xmax=xc+cube_side/2;
ymax=yc+cube_side/2;
zmax=zc+cube_side/2;
cube=X>=xmin & X<=xmax & Y>=ymin & Y<=ymax & Z>=zmin & Z<=zmax;
xcube=X(cube);
ycube=Y(cube);
zcube=Z(cube);

cube_rotated=zeros(3,x*numel(xcube));
% 
for k=1:x
    
R1=rotx(rand(1)*180)*roty(rand(1)*180)*rotz(rand(1)*180);
cube_rotated1=R1*([xcube,ycube,zcube]-repmat(mean([xcube,ycube,zcube])',[1,size(xcube,1)])')'+repmat(mean([xcube,ycube,zcube])',[1,size(xcube,1)]);
cube_rotated1=round(cube_rotated1);
cube_rotated1(3,:)=cube_rotated1(3,:)+M/((rand(1)-0.5)*60);
cube_rotated1(2,:)=cube_rotated1(2,:)+N/((rand(1)-0.5)*60);
cube_rotated1(1,:)=cube_rotated1(1,:)+N/((rand(1)-0.5)*60);
cube_rotated1=round(cube_rotated1);

cube_rotated(:,1+(k-1)*numel(xcube):k*numel(xcube))=cube_rotated1;

end

for i=1:size(cube_rotated,2)
    xx=cube_rotated(1,i);
    yy=cube_rotated(2,i);
    zz=cube_rotated(3,i);
    try
        xy_limit=xc+r_in(round(zz));
    catch
        xy_limit=xc+r_in(1);
    end
    if (xx>=1) && xx<=xy_limit &&  yy>=1 && yy<=xy_limit && zz>=1 && zz<=M
    phantomt(xx,yy,zz)=0.6;
    end
end


end

