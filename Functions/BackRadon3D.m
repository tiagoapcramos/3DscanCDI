function [ volume,proj_geom ] = BackRadon3D(frames,data)
% BACKRADON BackProjects the projections in 'frames' according to the
% projection parameters and volume geometry in data
%
% volume              - 3D volume / Tomographic Reconstruction
% data.parameters     - Matrix with projection parameters for all
%                       projections, in the form [Theta, u, v, alpha, beta]
%                       (:,1) - projection angles (degrees)
%                       (:,2) - horizontal granslations (Pixels)
%                       (:,3) - vertical translations (Pixels)
%                       (:,4) - alpha (Degrees)
%                       (:,5) - beta (Degrees)
% data.max_memory     - Maximum available memory in GPU
% frames              - 3D matrix with stack of frames
% proj_geom           - Structure array with 3D parallel beam geometry (see
%                       ASTRA's toolbox documentation)
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

% Input variables initialization
if ~isfield(data,'max_memory')
    max_memory=0.01;
else
    max_memory=data.max_memory;
end
if size(data.parameters,2)<5
    aux=zeros(size(data.parameters,1),5);
    aux(:,1:size(data.parameters,2))=data.parameters;
    data.parameters=aux;
end
sp=data.sp;
parameters=data.parameters;


% Create ASTRA volume geometry
% try
%     N=data.N;
%     M=data.M;
% catch
    N=ceil( max(parameters(:,2))-min(parameters(:,2)) +sp(2));
    M=ceil( max(parameters(:,3))-min(parameters(:,3)) +sp(1));
% end

vol_geom=astra_create_vol_geom(N,N,M);
volume=zeros(N,N,M);
%vol_geom=astra_create_vol_geom(N-sp(2),N-sp(2),M-sp(1));
%volume=zeros(N-sp(2),N-sp(2),M-sp(1));

% Create Astra's projector geometry
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

if isfield(data,'filter_type') && (strcmp(data.filter_type,'dc') || ...
        strcmp(data.filter_type,'f0'))
    for k=1:size(frames,3)
        frames(:,:,k)=frames(:,:,k)-squeeze(mean(mean(frames(:,:,k))));
    end
end

frames=permute(frames,[2,3,1]);

if isfield(data,'filter_type') && isfield(data,'filter_d')
    filt = DesignFilter(data.filter_type, data.sp(1), data.filter_d,0);
    H=repmat(filt,[1,size(frames,2),size(frames,1)]);
    frames(size(H,1),:,:)=0;%zero pad projection
    frames=fft(frames,[],1).*H;
    frames=ifft(frames,[],1);
    frames(data.sp+1:end,:,:)=[];
end

%Define forward and Backward operators
Atfun=@(fram) backproj(fram, proj_geom, vol_geom);
vol_memory=N*N*M*8/(1024^3);%number of voxels x 8bits x 2 (complex array)
req_memory=(sp(1)*sp(2)*size(parameters,1)+N*N*M)*8/(1024^3);
ir=isreal(frames);
if req_memory<max_memory
volume=Atfun(real(frames))+1i.*Atfun(imag(frames));
else
div=ceil((req_memory-vol_memory)/(max_memory-vol_memory));
if div<0,div=4;end
%div = ceil(req_memory/max_memory);
n_partial_proj=ceil(size(parameters,1)/div);
indexes=1:n_partial_proj;
	for k=1:div
		proj_geom_k=proj_geom;
		proj_geom_k.Vectors=proj_geom_k.Vectors(indexes,:);
		Atfun=@(fram) backproj(fram, proj_geom_k, vol_geom);
        try
            if ir==1
                volume=volume+Atfun((frames(:,indexes,:)));
            else
                volume=volume+Atfun(real(frames(:,indexes,:)))+1i.*Atfun(imag(frames(:,indexes,:)));
            end
        catch
            warning('This partition had no data to backproject...')
        end
		indexes=indexes+n_partial_proj;
		indexes(indexes>size(parameters,1))=[];
		astra_mex_data3d('clear');
	end
end

end


function [BP] = backproj(proj,proj_geom,vol_geom)
% Back Projection operator. Performs a Back-Projection operation from the
% projections proj, according to the projector geometry defined by 
% proj_geom. This operation is equivalent to A'*x
%
% This file is part of AutoTomoAlign, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) March-2017

projection_id=astra_mex_data3d('create','-sino',proj_geom,proj);
BP_id=astra_mex_data3d('create','-vol',vol_geom,0);
cfg=astra_struct('BP3D_CUDA');
cfg.ReconstructionDataId=BP_id;
cfg.ProjectionDataId=projection_id;
alg_id=astra_mex_algorithm('create',cfg);
astra_mex_algorithm('iterate',alg_id);
BP=astra_mex_data3d('get',BP_id);
astra_mex_algorithm('clear')
astra_mex_data3d('clear')

end

function filt = DesignFilter(filter_type, N, d,derivative)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


order = max(64,2^nextpow2(2*N));

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.


if derivative
    filt = 0*( 0:(order/2) )+1; %ORDER IS ALWAYS EVEN
else
    filt = 2*( 0:(order/2) )./order;
end
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter_type
    case 'none'
        filt=filt./( 2*( 0:(order/2) )./order);
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end

filt(w>pi*d) = 0;                      % Crop the frequency response
if derivative
    filt = [filt' ; -filt(end-1:-1:2)']/(1i*pi);    % Symmetry of the filter
else
    filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter IDEA FROM MANUEL GUIZAR Sicairos
end
end