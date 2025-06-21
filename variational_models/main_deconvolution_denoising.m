%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Image Deconvolution and Denoising Analysis                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Imports
clear all, close all, clc

scriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptPath);
commonsPath = fullfile(projectRoot, 'commons');

addpath(commonsPath);

%% Read Image

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'images');
im = imread(fullfile(img_folder, 'panda.png'));
im = im2double(im); % im2double transforms the data into double type and in the range [0,1]


%% Blurring with Noise
% Gaussian Noise
Noise_STD = 0;
noise_im = Noise_STD*randn(size(im));

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
im_blur_noise = R(im) + noise_im;

set(0, 'DefaultFigureColor', 'w')

% Visualizacion
figure('Color', 'w') 
subplot(141), imagesc(im), title('Original image')
axis off
axis image
subplot(142), MAP=colormap; imagesc(kernel),  title('kernel')
axis off
axis image
subplot(143), MAP=colormap; imagesc(noise_im), title('Image Noise')
axis off
axis image
subplot(144), imagesc(im_blur_noise), title('Blurred + Noise image')
axis off
axis image


%% DECONVOLUTION & DENOISING ALGORITHM

% Algorithm input parameters
varin.lambda1   = 1;              % Fidelity Hyperparam
varin.lambda2   = 1;              % Fidelity Hyperparam
varin.Nit       = 10000;          % Number of Iterations
varin.dt        = 1*1e-2;         % Step size 
varin.epsilon   = 1.e-6;          % Epsilon
varin.hasNoise  = false;          % Image has Noise

varin.f         = im_blur_noise;  % Blurred image
varin.im_org    = im;             % Original Image
varin.kernel    = kernel;         % Blurred Kernel
varin.p         = 2;              % p Laplacian
varin.Verbose   = 2;              % Verbose

% Algorithm execution
[varout] = Deconvolution_Denoising_model(varin);

% Algortihm output
u_p = varout.u;

if varin.Verbose ~= 0
    psnr_p = varout.psnr;
    en_p   = varout.en;
    pr_p   = varout.pr;
    fi_p   = varout.fi;
    ssim_p = varout.ssim;
end

%% Show Deconvolution Model Final Results
[ssimval,ssimmap] = ssim(u_deconv,im);

figure
subplot(141), imagesc(u_p), title(['Reconstructed Image PSNR: ',num2str(PSNR(im,u_deconv)),' db'])
axis off
axis image
subplot(142), plot(en_p,'r','LineWidth', 2), hold on,
              plot(fi_p,'b', 'LineWidth', 2),
              plot(pr_p,'g','LineWidth', 2)
              legend('Total Energy','Fidelity','Prior'), grid on
subplot(143), plot(psnr_p,'c','LineWidth', 2), 
              legend('PSNR (db)'), grid on
subplot(143), plot(ssim_p,'c','LineWidth', 2), 
              legend('SSIM'), grid on
subplot(144), imagesc(ssimmap), title("Reconstructed Image SSIM: "+ssimval)
axis off
axis image

%% Original Image vs Reconstructed Comparisson

figure,
subplot(131), imagesc(im),       title('Original Image')
axis off
axis image
subplot(132), imagesc(im_blur_noise),  title('Blurred Image')
axis off
axis image
subplot(133), imagesc(u_p),      title('Restored Image')
axis off
axis image


%% Matlab Deconvolution Models: Lucy-Richardson deconvolution | Deconvreg

% Lucy-Richardson Model
lucy = deconvlucy(im_blur_noise, kernel, 100);
[ssim_lucy,ssimmap_lucy] = ssim(lucy,im);

% DeconvReg model
[u_deconvreg, LAGRA] = deconvreg(im_blur_noise, kernel, [], 1);
[ssim_deconvreg,ssimmap_deconvreg] = ssim(u_deconvreg,im);

% Visalization
figure,
subplot(221), imagesc(u_deconvreg),title(['Deconvreg-L1 PSNR: ',num2str(PSNR(im,u_deconvreg)),' db'])
axis off
axis image
subplot(224), imagesc(ssimmap_deconvreg), title("Deconvreg-L1 SSIM: "+ssim_lucy)
axis off
axis image
subplot(223), imagesc(lucy),    title(['Lucy-Richardson PSNR: ',num2str(PSNR(im,lucy)),' db'])
axis off
axis image
subplot(224), imagesc(ssimmap_lucy), title("Lucy-Richardson SSIM: "+ssim_lucy)
axis off
axis image


%% Final Comparisson with Matlab Methods

% Visalization
figure, 
subplot(141), imagesc(im),       title('Original image')
axis off
axis image
subplot(142), imagesc(u_deconvreg), title(['Deconvreg-L1: ', num2str(PSNR(im,u_deconvreg)),' db'])
axis off
axis image
subplot(143), imagesc(luc1),     title(['Lucy-Richardson: ', num2str(PSNR(im,lucy)),' db'])
axis off
axis image
subplot(144), imagesc(u_p), title(['p=',num2str(varin.p),': ', num2str(PSNR(im,u_p)),' db'])
axis off
axis image