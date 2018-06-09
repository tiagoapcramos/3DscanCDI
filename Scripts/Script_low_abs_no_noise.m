close all
clear 
clc
%% %%%       MAIN SCRIPT FOR 3D PTYCHOGRAPHIC RECONSTRUCTION      %%% %%
% This script is used to generate simulated ptychography data and solve a
% combined 3D phase-retrieval and tomographic reconstruction of far-field
% diffraction patterns. This non-linear inverse problem is solved by means
% of the Levenberg-Marquadt algorithm. 
% The script is organized in the following sections:
%   1   - Add relevant paths (functions and data) to your workspace
%   2   - Generate Experimental Setup, Phantom, probe function and
%         ptychography data
%   3   - Define Solver parameters
%   4   - Solve the Inverse Problem
%   5   - Monitor, image results and save
%
% Tiago Ramos (tiagoj@dtu.dk) 03/01/2018
%% 1 - ADD RELEVANT PAHTS
astra_path={'C:\ASTRA\astra-1.6\mex','C:\ASTRA\astra-1.6\tools'};
% astra_path={'/zhome/a5/2/93279/astra-1.6/matlab/mex',...
%     '//zhome/a5/2/93279/astra-1.6/matlab/tools'};
tools_path={'../Functions','../External','../Data'};
addpath(astra_path{:},tools_path{:});
%% 2 - Experimental Setup, and Data Generation
%Size of the test object
N = 300;
%Number of diffraction patterns to simulate
n_dp = 4000;
%Add Poisson noise
add_noise = false;
%High absorbing material
high_absorption = false;
%Generates Simulated Phantom and experimental parameters
[volume, data]=CreateTestProblem(N,n_dp);
%Normalize volume to chosen values
volume = real(volume)*1e-5+1i.*imag(volume).*1e-7;
if high_absorption
    volume = real(volume)+1i.*imag(volume).*1e2;
end
% Compute Diffraction patterns. Round to integer numbers and add noise
data.Imeas = floor(ForwardF(volume,data));
if add_noise
    data.Imeas = imnoise(data.Imeas *1e-12,'poisson')*1e12;
end
%% 3. Define solver parameters
options.mu0 = NaN;
options.n_iter = 49;
options.n_iter_cgm = 24;
options.noise_model = 'poisson';
options.monitor = true;

%% 4. Solve Inverse Problem
[recon, output]=Reconstruct_LMA(data,options);
%% 5. Monitor results and save
figure('color','w')
SliceViewer(volume)
colormap(flipud(bone));
figure('color','w')
SliceViewer(recon)
colormap(flipud(bone))

s = whos;
for i = 1:length(s)
      if strcmp(s(i).class,'double')
          name = s(i).name;
          assignin('base', name, single(evalin('base', name)));
      end
end

save('results_ML.mat','-v7.3');