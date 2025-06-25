%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Convolution Analytic Reconstruction                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Imports
clear all, close all, clc

basePath = fileparts(pwd);
commonsPath = fullfile(basePath, 'commons');
addpath(commonsPath);

%% Read Image
clear all, close all, clc

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'imagenes');
im = imread(fullfile(img_folder, 'panda.png'));
im = im2double(im); % im2double transforms the data into double type and in the range [0,1]


%% Blurring with Noise
% Gaussian Noise
Noise_STD = 0;

% Load Gaussian Kernel
kernel_size = 21;
sigma = 5;
kernel_gaussian = fspecial('gaussian', [kernel_size, kernel_size], sigma);

% Load Motion Kernels
load motion_kernels.mat;

kernel = kernel7; %kernel_gaussian

% Kernel and Fourier Transform Kernel
dim = size(im);
kernel_Fourier = psf2otf(kernel,[dim(1),dim(2)]);

% Image Transformation Kernel + Noise
R =@(x) real(ifft2(kernel_Fourier.*fft2(x)));
im_blur_noise = R(im) + Noise_STD*randn(dim);


%% IMAGE RECONSTRUCTION (ANALITYC WAY 1)
% Based on the image convolution formula
% f = k * u -> u = F^-1(F(f)/F(k))

u_analytic = real(ifft2(fft2(im_blur_noise)./kernel_Fourier));

psnr_u = PSNR(u_analytic,im);

% Visualizacion
set(0, 'DefaultFigureColor', 'w')

figure('Color', 'w') 
subplot(131), imagesc(im), title(['Original image'])
axis off
axis image
subplot(132), imagesc(im_blur_noise), title(['Blurred image'])
axis image
axis off
subplot(133), imagesc(u_analytic),  title(['Analytic Reconstructed Image PSNR: ',num2str(psnr_u)])
axis image
axis off


%% IMAGE RECONSTRUCTION (ANALITYC WAY 2)
% Based on the Euler-Lagrange Equations
% k^T * ( k * u - f ) = 0 -> u = F^-1( ( F(k^T) · F(f) ) / | F(k) |^2 )

u_reconstructed = real(ifft2(conj(kernel_F).*fft2(im_blur_noise)./(abs(kernel_Fourier).^2)));

psnr_u = PSNR(u_reconstructed,im);

% Visualizacion
set(0, 'DefaultFigureColor', 'w')

figure('Color', 'w') 
subplot(131), imagesc(im), title(['Original image'])
axis off
axis image
subplot(132), imagesc(im_blur_noise), title(['Blurred image'])
axis image
axis off
subplot(133), imagesc(u_reconstructed),  title(['Analytic Reconstructed Image PSNR: ',num2str(psnr_u)])
axis image
axis off