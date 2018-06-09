function [volume,data ] = CreateTestProblem(N,n_dp)
%CREATETESTPROBLEM - generates 3D phantom and probe function for the 3D
%ptychography simulations.
%Input:
% N     -   Size of the simulated object: cube of dimensions [NxNxN]
% n_dp  -   number of diffraction patterns to simulate
%
% This file is part of 3DPtychoTomo, which is released under the
% BSD 3-Clause License. Please see LICENSE.txt
% Tiago Ramos (tiagoj@dtu.dk) May-2018

data.N = N;
data.M = N;
%probe size in pixels
data.sp = floor([N,N]*1/3);
%X-ray wavelength [m]
data.lambda = 1e-10;
%Sample-detector distance [m]
data.z = 5;
%Detector pixel size [m]
detect_pix_size = 172e-6;
%Reconstructed pixel size
data.pixsize = data.z*data.lambda/(data.sp(1)*detect_pix_size);
%Wave number
data.k = (2*pi*data.pixsize)/data.lambda;
%Number of projection angles. Here i fix the 16, meaning that for each
% projection angle we will have 16 diffraction patterns randomly oriented
data.n_proj = ceil(n_dp/16);
%Maximum avaliable GPU memory [GB]
data.max_memory = 1;
%Assign amplitude of random orientation parameters
uncert_options.theta_u = 0;
uncert_options.u_u = data.N/2;
uncert_options.v_u = data.M/2;
[data.parameters, data.proj_geom]=CreateRandomParameters(data.n_proj*16,...
    uncert_options,data.sp);
%Create Probe function from Complex Gaussian function
data.probe = GaussianProbe(data.sp,[],[],1);
%Mask Probe Function and normalize its intensities based on flux and
%exposure time
data.probe = data.probe.*Mask_Circle(round(0.9*data.sp(1)),data.sp(1));
exposure_time = 0.1;%[s]
total_flux = 1e8; %[photons/s]
data.probe = data.probe * sqrt(total_flux*exposure_time/...
    sum(abs(data.probe(:)).^2));

%Generate volume from [0 1]
%The real and imaginary parts of the generated phantom should be properly
%scaled afterwards
volume = newphantom2(data.N,data.M,100);
volume = volume/max(volume(:));
volume2 = phantom3d_better(data.N);
volume = -volume + 1i.*volume2;

end

