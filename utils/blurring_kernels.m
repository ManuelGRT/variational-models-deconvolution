%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Types of Blurring Kernels                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read Image
clear all, close all, clc

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'images');
im = imread(fullfile(img_folder, 'panda.png'));
im = im2double(im); % im2double transforms the data into double type and in the range [0,1]

%% Gaussian Blurring with Noise
% Gaussian Noise
Noise_STD = 0;

% Load Gaussian Kernel
kernel_size = 21;
sigma = 5;
kernel_gaussian = fspecial('gaussian', [kernel_size, kernel_size], sigma);
kernel = kernel_gaussian;

% Kernel and Fourier Transform Kernel
dim = size(im);
kernel_Fourier = psf2otf(kernel,[dim(1),dim(2)]);

% Image Transformation Kernel + Noise
R =@(x) real(ifft2(kernel_Fourier.*fft2(x)));
im_blur_noise = R(im) + Noise_STD*randn(dim);

set(0, 'DefaultFigureColor', 'w')

% Visualizacion
figure('Color', 'w') 
subplot(131), imagesc(im), title(['Original image'])
axis off
axis image
subplot(132), MAP=colormap; imagesc(kernel),  title('Gaussian kernel')
axis off
axis image
subplot(133), imagesc(im_blur_noise), title(['Blurred image'])
axis off
axis image


%% Motion Blurring with Noise
% Gaussian Noise
Noise_STD = 0;

% Load Motion Kernels
load kernels.mat;
kernel = kernel7;

% Kernel and Fourier Transform Kernel
dim = size(im);
kernel_Fourier = psf2otf(kernel,[dim(1),dim(2)]);

% Image Transformation Kernel + Noise
R =@(x) real(ifft2(kernel_Fourier.*fft2(x)));
im_blur_noise = R(im) + Noise_STD*randn(dim);

set(0, 'DefaultFigureColor', 'w')

% Visualizacion
figure('Color', 'w') 
subplot(131), imagesc(im), title(['Original image'])
axis off
axis image
subplot(132), MAP=colormap; imagesc(kernel),  title('Motion kernel')
axis off
axis image
subplot(133), imagesc(im_blur_noise), title(['Blurred image'])
axis off
axis image
