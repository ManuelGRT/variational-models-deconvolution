%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Types of Image Convolution                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read Image
clear all, close all, clc

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'images');
im = imread(fullfile(img_folder, 'panda.png'));
im = im2double(im); % im2double transforms the data into double type and in the range [0,1]


%% Image Convolutions in MATLAB

% Load Gaussian Kernel
kernel_size = 21;
sigma = 5;
kernel_gaussian = fspecial('gaussian', [kernel_size, kernel_size], sigma);

% Load Motion Kernels
load motion_kernels.mat;

kernel = kernel7; %kernel_gaussian

% Empty structures for different convolutions
kernel_size = [size(kernel), 1];
im_blur_full  = zeros(size(im) + kernel_size - 1);
im_blur_same  = zeros(size(im));
im_blur_valid = zeros(size(im,1)-size(kernel,1)+1, size(im,2)-size(kernel,2)+1, 3);

% Apply convolution by channel
for c = 1:3
    u = im(:,:,c);               % Canal R, G, B

    % Full padding
    im_blur_full(:,:,c) = conv2(u, kernel, 'full');

    % Same padding (Same size)
    im_blur_same(:,:,c) = conv2(u, kernel, 'same');

    % Valid padding
    im_blur_valid(:,:,c) = conv2(u, kernel, 'valid');
end


% Kernel and Convolutional Image with Fourier Transform
dim = size(im);
kernel_Fourier = psf2otf(kernel,[dim(1),dim(2)]);
R =@(x) real(ifft2(kernel_Fourier.*fft2(x)));
im_blur_Fourier = (R(im));


% Visualización
set(0, 'DefaultFigureColor', 'w')
figure('Color', 'w') 

subplot(1,5,1), imagesc(im),  title('Original Image')
axis off
axis image
subplot(1,5,2), imagesc(im_blur_full), title('Blurring: conv2(..., ''full'')')
axis off
axis image
subplot(1,5,3), imagesc(im_blur_same), title('Blurring: conv2(..., ''same'')')
axis off
axis image
subplot(1,5,4), imagesc(im_blur_valid), title('Blurring: conv2(..., ''valid'')')
axis off
axis image
subplot(1,5,5), imagesc(im_blur_Fourier), title('Blurring: Fourier Transform')
axis off
axis image
