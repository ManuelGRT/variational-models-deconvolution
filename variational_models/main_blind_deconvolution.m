%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Image Blind Deconvolution Analysis                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Imports
clear all, close all, clc

basePath = fileparts(pwd);
commonsPath = fullfile(basePath, 'commons');
addpath(commonsPath);

%% Read Image

basePath = fileparts(pwd);
img_folder = fullfile(basePath, 'images');
im = imread(fullfile(img_folder, 'flor.jpeg'));
im = im2double(im); % im2double transforms the data into double type and in the range [0,1]


%% Image Blurring

% Load Gaussian Kernel
kernel_size = 21;
sigma = 5;
kernel_gaussian = fspecial('gaussian', [kernel_size, kernel_size], sigma);

% Load Motion Kernels
load motion_kernels.mat;

kernel = kernel8; %kernel_gaussian

% Image Convolution
im_blur_f = convn(im,kernel,'valid');
im_blur_u = convn(im,kernel,'same');

set(0, 'DefaultFigureColor', 'w')

% Visualizacion
figure('Color', 'w') 
subplot(131), imagesc(im), title(['Original Image'])
axis off
axis image
subplot(132), MAP=colormap; imagesc(kernel),  title(['Kernel'])
axis off
axis image
subplot(133), imagesc(im_blur_u), title(['Blurred Image'])
axis off
axis image


%% DECONVOLUTION ALGORITHM

% Algorithm input parameters
varin.Nit           = 10000;            % Number of iterations
varin.epsilon_stop  = 1.e-6;            % Epsilon stop criteria
varin.dt_u          = 1*1e-3;           % Step size u
varin.dt_k          = 1*1e-3;           % Step size Kernel
varin.lambda        = 2.5;                % Regularization Hyperparam
varin.p             = 1;                % p Laplacian

varin.f             = im_blur_f;        % Blurred Image
varin.u             = im_blur_u;        % Initial u
varin.im_org        = im;               % Original Image
varin.kernel        = kernel;           % Blurred Kernel
varin.Verbose       = 2;                % Verbose

% Algorithm execution
[varout] = Blind_Deconvolution_model(varin);

% Algortihm output
u_blind    = varout.u;
k_blind    = varout.k;

if varin.Verbose ~= 0
    psnr_blind = varout.psnr;
    en_blind   = varout.en;
    pr_blind   = varout.pr;
    fi_blind   = varout.fi;
    ssim_blind = varout.ssim_u;
end

%% Show Blind Deconvolution Model Final Results
[ssimval,ssimmap] = ssim(u_blind,im);

figure
subplot(161), imagesc(u_blind), title(['Reconstructed Image PSNR: ',num2str(PSNR(im,u_blind)),' db'])
axis off
axis image
subplot(162), plot(en_blind,'r','LineWidth', 2), hold on,
              plot(fi_blind,'b', 'LineWidth', 2),
              plot(pr_blind,'g','LineWidth', 2)
              legend('Total Energy','Fidelity','Prior'), grid on
subplot(163), plot(psnr_blind,'c','LineWidth', 2),
              legend('PSNR (db)'), grid on
subplot(164), plot(ssim_blind,'c','LineWidth', 2), 
              legend('SSIM'), grid on
subplot(165), imagesc(ssimmap), title("Reconstructed Image SSIM: "+ssimval)
axis off
axis image
subplot(166), imagesc(k_blind), title('Reconstructed kernel ')
axis off
axis image

%% Original Image vs Reconstructed Comparisson

figure,
subplot(131), imagesc(im),       title('Original Image')
axis off
axis image
subplot(132), imagesc(im_blur_u),  title('Blurred Image')
axis off
axis image
subplot(133), imagesc(u_blind), title('Restored Image')
axis off
axis image


%% Original Kernel vs Reconstructed Comparisson

figure,
subplot(121), MAP=colormap; imagesc(kernel),  title('Original kernel')
axis off
axis image
subplot(122), MAP=colormap; imagesc(k_blind),  title('Reconstructed kernel')
axis off
axis image


%% Matlab Blind Deconvolution Model: DeconvBlind

% DeconvBlind
init_kernel =  fspecial('gaussian', [21, 21], 0.001);  % Initial Kernel
[u_deconvblind, k_deconvblind] = deconvblind(im_blur_u, init_kernel, 10);
[ssim_deconvblind,ssimmap_deconvblind] = ssim(u_deconvblind,im);


% Visalization
figure,
subplot(131), imagesc(u_deconvblind),title(['DeconvBlind PSNR: ',num2str(PSNR(im,u_deconvblind)),' db'])
axis off
axis image
subplot(132), imagesc(ssimmap_deconvblind), title("DeconvBlind SSIM: "+ssim_deconvblind)
axis off
axis image
subplot(133), imagesc(k_deconvblind),    title('DeconvBlind Reconstructed Kernel')
axis off
axis image


%% Final Image Comparisson with Matlab Methods

% Visalization
figure, 
subplot(131), imagesc(im),       title('Original image')
axis off
axis image
subplot(132), imagesc(u_deconvblind), title(['DeconvBlind: ', num2str(PSNR(im,u_deconvblind)),' db'])
axis off
axis image
subplot(133), imagesc(u_blind), title(['Blind Deconvolution: ', num2str(PSNR(im,u_blind)),' db'])
axis off
axis image


%% Final Kernel Comparisson with Matlab Methods

% Visalization
figure, 
subplot(131), imagesc(kernel),          title('Original kernel')
axis off
axis image
subplot(132), imagesc(k_deconvblind),   title('DeconvBlind Reconstructed Kernel')
axis off
axis image
subplot(133), imagesc(k_blind),         title('Blind Deconvolution Reconstructed Kernel')
axis off
axis image