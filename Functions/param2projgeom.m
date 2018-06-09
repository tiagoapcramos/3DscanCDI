function [ proj_geom ] = param2projgeom( parameters,sp )
%PARAM2PROJGEOM creates an ASTRA projector geometry (proj_geom) from an
%array of projection parameters (parameters)
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018


%Define ASTRA projection geometry
proj_geom=astra_create_proj_geom('parallel3d',1,1,sp(1),...
    sp(2),parameters(:,1)*pi/180);
proj_geom=astra_geom_2vec(proj_geom);

proj_geom.DetectorRowCount=sp(1);
proj_geom.DetectorColCount=sp(2);
if size(parameters,2)<=3 || sum(sum(parameters(4:end,1:end)))==0
	proj_geom.Vectors(:,1)=sind(parameters(:,1));                       %rayX
	proj_geom.Vectors(:,2)=-cosd(parameters(:,1));                      %rayY
	proj_geom.Vectors(:,3)=0;                                           %rayZ
	proj_geom.Vectors(:,4)=(parameters(:,2)).*cosd(parameters(:,1));    %dX
	proj_geom.Vectors(:,5)=(parameters(:,2)).*sind(parameters(:,1));    %dY
	proj_geom.Vectors(:,6)=(parameters(:,3));                           %dZ
	proj_geom.Vectors(:,7)=cosd(parameters(:,1));                       %uX
	proj_geom.Vectors(:,8)=sind(parameters(:,1));                       %uY
	proj_geom.Vectors(:,9)=0;                                           %uZ
	proj_geom.Vectors(:,10)=0;                                          %vX
	proj_geom.Vectors(:,11)=0;                                          %vY
	proj_geom.Vectors(:,12)=1;                                          %vZ
else
	% Define 3D Transformation/rotation matrix
	for i=1:size(parameters,1)
		R=rotz(parameters(i,1))*rotx(parameters(i,4))*roty(parameters(i,5));
		u(i,:)=(R*[1 0 0]')';                     
		v(i,:)=(R*[0 0 1]')';   
	end                    
	proj_geom.Vectors(:,7:12)=[u v];   
	proj_geom.Vectors(:,4:5)=proj_geom.Vectors(:,7:8).*...
		[parameters(:,2) parameters(:,2)]; 
	proj_geom.Vectors(:,6)=proj_geom.Vectors(:,12).*parameters(:,3);
	proj_geom.Vectors(:,1:3)=cross(proj_geom.Vectors(:,7:9)',...
		proj_geom.Vectors(:,10:12)')';  

end

end

