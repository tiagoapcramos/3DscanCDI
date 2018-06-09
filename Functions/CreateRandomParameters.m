function [ parameters,proj_geom] = CreateRandomParameters(num_proj,...
    uncert_options,sp)
% CREATERANDOMPARAMETERS Creates an array with random projection parameters
% define by the uncertanty options

%   num_proj -  Number of projection angles
%   uncert_options - Struct variable to define errors in the projection
%                   parameters
%        - theta_u - amplitude of random distributed errors in theta
%        - u_u - amplitude of random distributed errors in u
%        - v_u - amplitude of random distributed errors in v
%        - alpha_u - amplitude of random distributed errors in alpha
%        - beta_u - amplitude of random distributed errors in beta
%        - rngs   - integer vector of length 5. Define 'seeds' for
%                   random number generation. Useful for data repeatability
%        - impulse - Bolean variable to include impulse noise in the
%                       translation parameters
%        - noise_perc - density of impulse noise [0 1]
%        - noise_strength - Amplitude of impulse noise
%
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

if nargin<2||isempty(uncert_options),uncert_options=struct;end
if nargin<3||isempty(sp),sp=[1,1];end

% Check for default values and missing fields in uncert_options
uncert_options=checkdefaultvalues(uncert_options);

angles_measured=linspace(0,180-180/num_proj,num_proj);

%Define theta errors
rng(uncert_options.rngs(1));
theta_errors = 2*uncert_options.theta_u*(rand(num_proj,1) - 0.5);
angles_true = angles_measured + theta_errors';

%Define u errors
rng(uncert_options.rngs(2));
u_errors=2*uncert_options.u_u*(rand(num_proj,1)-0.5);
if uncert_options.impulse==true
    rng(uncert_options.rngs(1));
    y_n1=imnoise(zeros(num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    rng(uncert_options.rngs(2));
    y_n2=-imnoise(zeros(num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    u_errors=u_errors+uncert_options.u_u*(y_n1+y_n2)*...
        uncert_options.noise_strength;
end

%Define v errors
rng(uncert_options.rngs(3));
v_errors=2*uncert_options.v_u*(rand(num_proj,1)-0.5);
if uncert_options.impulse==true
    rng(uncert_options.rngs(3));
    y_n1=imnoise(zeros(num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    rng(uncert_options.rngs(4));
    y_n2=-imnoise(zeros(num_proj,1),'salt & pepper',...
        uncert_options.noise_perc/100);
    v_errors=v_errors+uncert_options.v_u*(y_n1+y_n2)*...
        uncert_options.noise_strength;
end

%Define alpha errors
rng(uncert_options.rngs(4));                                                                      
alpha_errors=2*uncert_options.alpha_u*(rand(num_proj,1)-0.5);

%Define beta errors                                                                       
rng(uncert_options.rngs(5));
beta_errors=2*uncert_options.beta_u*(rand(num_proj,1)-0.5);

%Define Projection parameters array
parameters=[angles_true', u_errors, v_errors, alpha_errors, beta_errors];

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

function [uncert_options] = checkdefaultvalues(uncert_options)
% Generates struct variable with default values for uncertanties options in
% case some or all are missing
if ~isfield(uncert_options,'theta_u'), uncert_options.theta_u=0;end
if ~isfield(uncert_options,'u_u'), uncert_options.u_u=0;end
if ~isfield(uncert_options,'v_u'), uncert_options.v_u=0;end
if ~isfield(uncert_options,'alpha_u'), uncert_options.alpha_u=0;end
if ~isfield(uncert_options,'beta_u'), uncert_options.beta_u=0;end
if ~isfield(uncert_options,'rngs'), uncert_options.rngs=[6 3 4 1 2];end
if ~isfield(uncert_options,'impulse'), uncert_options.impulse=false;end
if ~isfield(uncert_options,'noise_perc'), uncert_options.noise_perc=0.1;end
if ~isfield(uncert_options,'noise_strength')
    uncert_options.noise_strength=.5;
end
end
