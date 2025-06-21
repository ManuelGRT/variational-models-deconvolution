%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Fourier Transform of an Image                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read Image
clear all, close all, clc

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'imagenes');
im = imread(fullfile(img_folder, 'iceberg_gris.jpg'));
im_gray = im2double(rgb2gray(im)); 
im_gray = im2double(im_gray); % im2double transforms the data into double type and in the range [0,1]

%% Fourier Transform of an Image

im_Fourier = fft2(im_gray);           % FFT
im_Fourier_centered = fftshift(im_Fourier);       % Centrar el cero
im_Fourier_Final = log(1 + abs(im_Fourier_centered));   % Escala logarítmica para mejor visualización

%% Show Figures

set(0, 'DefaultFigureColor', 'w')

figure
subplot(121), imshow(im), title(['Original image'])
axis off
subplot(122), imshow(im_Fourier_Final,[]), title(['Fourier Transform image'])
axis off