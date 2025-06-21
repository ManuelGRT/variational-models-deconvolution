%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Image Deconvolution Analysis                       %
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


%% Image Blurring

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
im_blur = R(im);

set(0, 'DefaultFigureColor', 'w')

% Visualizacion
figure('Color', 'w') 
subplot(131), imagesc(im), title(['Original Image'])
axis off
axis image
subplot(132), MAP=colormap; imagesc(kernel),  title(['Kernel'])
axis off
axis image
subplot(133), imagesc(im_blur), title(['Blurred Image'])
axis off
axis image


%% DECONVOLUTION ALGORITHM

% Algorithm input parameters
varin.Nit           = 10000;            % Number of iterations
varin.epsilon_stop  = 1.e-6;            % Epsilon stop criteria
varin.dt            = 1*1e-2;           % Step size

varin.f             = im_blur;          % Blurred Image
varin.im_org        = im;               % Original Image
varin.kernel        = kernel;           % Blurred Kernel
varin.Verbose       = 2;                % Verbose

% Algorithm execution
[varout] = Deconvolution_model(varin);

% Algortihm output
u_deconv    = varout.u;

if varin.Verbose ~= 0
    psnr_deconv = varout.psnr;
    en_deconv   = varout.en;
    pr_deconv   = varout.pr;
    fi_deconv   = varout.fi;
    ssim_deconv = varout.ssim;
end

%% Show Deconvolution Model Final Results
[ssimval,ssimmap] = ssim(u_deconv,im);

figure
subplot(141), imagesc(u_deconv), title(['Reconstructed Image PSNR: ',num2str(PSNR(im,u_deconv)),' db',',  SSIM: ',num2str(ssimval)])
axis off
axis image
subplot(142), plot(en_deconv,'r','LineWidth', 2), hold on,
              plot(fi_deconv,'b', 'LineWidth', 2),
              plot(pr_deconv,'g','LineWidth', 2)
              legend('Total Energy','Fidelity','Prior'), grid on
subplot(143), plot(psnr_deconv,'c','LineWidth', 2),
              legend('PSNR (db)'), grid on
subplot(143), plot(ssim_deconv,'c','LineWidth', 2), 
              legend('SSIM'), grid on
subplot(144), imagesc(ssimmap), title("Reconstructed Image SSIM: "+ssimval)
axis off
axis image

%% Original Image vs Reconstructed Comparisson

figure,
subplot(131), imagesc(im),       title('Original Image')
axis off
axis image
subplot(132), imagesc(im_blur),  title('Blurred Image')
axis off
axis image
subplot(133), imagesc(u_deconv), title('Restored Image')
axis off
axis image


%% Matlab Deconvolution Models: Lucy-Richardson deconvolution | Deconvreg

% Lucy-Richardson Model
lucy = deconvlucy(im_blur, kernel, 100);
[ssim_lucy,ssimmap_lucy] = ssim(lucy,im);

% DeconvReg model
[u_deconvreg, LAGRA] = deconvreg(im_blur, kernel, [], 1);
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
subplot(144), imagesc(u_deconv), title(['Deconvolution: ', num2str(PSNR(im,u_deconv)),' db'])
axis off
axis image