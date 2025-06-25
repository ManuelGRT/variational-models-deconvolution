%% --------------------------------------------------------------------- %%
%      Blind Deconvolution: Perrone, Favaro, E. Schiavi, M. Ramirez       %
%%-----------------------------------------------------------------------%%
function [varout]= Blind_Deconvolution_model(varin)
%-------------------------------------------------------------------------%
%                                 Inputs                                  %
%-------------------------------------------------------------------------%
im_org   = varin.im_org;
f        = varin.f;
u        = varin.u;
kernel   = varin.kernel;

p        = varin.p;
Nit      = varin.Nit;
dt_u     = varin.dt_u;
dt_k     = varin.dt_k;
lambda   = varin.lambda;
epsilon  = varin.epsilon_stop;
Verbose  = varin.Verbose;

switch Verbose
    case 0

    case 1

    case 2
        figure,
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15, 0.3, 0.7, 0.4]);
end

%-------------------------------------------------------------------------%
%                               Algorithm                                 %
%-------------------------------------------------------------------------%

%--------------------- Kernel & Image dimensions -------------------------%
[MU, NU, C] = size(im_org);
[M, N, C] = size(f);
[MK,NK]   = size(kernel);

%------------------------ Convolution Operations -------------------------%
Conv_u  =@(u,k) conv2(u, k, 'valid');
Convt_u =@(u,k) conv2(u, rot90(k,2), 'full');

Convt_k =@(u,f) conv2(rot90(u,2), f, 'valid');

% Kernel and Fourier Transform Kernel
dim = size(im_org);
kernel_Fourier = psf2otf(kernel,[dim(1),dim(2)]);

% Image Transformation Kernel + Noise
R =@(x) real(ifft2(kernel_Fourier.*fft2(x)));
f_same = u;
%----------------------- Variables initialization ------------------------%
k = ones(MK,NK)/MK/NK;

iter = 1;
stop = false;

while (iter < Nit)  && (~stop)
%---------------------------- Visualization ------------------------------%
    switch Verbose
        case 0

            disp(['Iter: ',num2str(iter)])
            
        case 1
            if ~mod(iter,100)
                imshow(u),pause(0.01)
            end
            [en(iter),pr(iter),fi(iter)] = pEnergy_R(u,f_same,lambda,p,R);
            psnr(iter) = PSNR(u,im_org);
            [ssim_u(iter),ssimmap] = ssim(u,im_org);
            disp(['Iter: ',num2str(iter),' PSNR: ', num2str(psnr(iter)),' SSIM: ', num2str(ssim_u(iter)),' Total Energy: ', num2str(en(iter)), ' Prior: ', num2str(pr(iter)), ' Fidelity: ',num2str(fi(iter)) ])
        case 2
            [en(iter),pr(iter),fi(iter)] = pEnergy_R(u,f_same,lambda,p,R);
            psnr(iter) = PSNR(u,im_org);
            [ssim_u(iter),ssimmap] = ssim(u,im_org);
            subplot(151), imshow(u), title(['PSNR: ',num2str(psnr(iter))])
            subplot(152), plot(en,'r','LineWidth', 2), hold on,
                          plot(fi,'b','LineWidth', 2),
                          plot(pr,'g','LineWidth', 2)
                          legend('Total Energy','Fidelity','Prior'), grid on
            subplot(153), plot(psnr,'c','LineWidth', 2), 
                          legend('PSNR (db)'), grid on
            subplot(154), imshow(ssimmap,[]), title("SSIM: "+ssim_u(iter))
            subplot(155), imagesc(k), title(['Reconstructed kernel'])
            axis image
            axis off
            pause(0)
    end

%----------------------------- Algorithm ---------------------------------%
    
    % Image Gradient
    ux = parcial_x(u);
    uy = parcial_y(u);

    % Modulus of the Image Gradient
    mod_grad = sqrt(ux.^2 + uy.^2 + epsilon^2);
    b        = mod_grad.^(p-2);

    % Image Laplacian   
    plap = div_x(b.*ux) + div_y(b.*uy);

    % Image update
    u_conv = zeros(MU,NU,C);
    for c=1:C
        u_conv(:,:,c) = Convt_u(Conv_u(u(:,:,c), k) - f(:,:,c), k);
    end
    ulap = u_conv - lambda*plap;
    dt_u = 5e-3*max(u(:))/max(1e-31, max(max(abs(ulap(:)))));
    u2   = u - dt_u*(ulap);

    % Kernel update
    k_conv = zeros(MK,NK);
    for c=1:C
        k_conv = k_conv + Convt_k( u(:,:,c), Conv_u(u(:,:,c), k) - f(:,:,c) );
    end
    dt_k = 1e-3*max(k(:))/max(1e-31, max(max(abs(k_conv))));
    k = k - dt_k*k_conv;
    
    % Kernel projection
    k = k.*(k>0);
    k = k/sum(k(:));

    % Stop criteria
    norm_u = norm(u2(:) - u(:));
    disp(['Stop Criteria: ',num2str(norm_u)])
    stop = norm_u < epsilon;

    % Variables update
    u = u2;
    iter = iter +1;
    lambda = max(1e-4,lambda * 0.999);

end
%--------------------------- Output Variables ----------------------------%
if Verbose == 0
    varout.u        = u;
    varout.k        = k;
else
    varout.u        = u;
    varout.k        = k;
    varout.en       = en;
    varout.pr       = pr;
    varout.fi       = fi;
    varout.psnr     = psnr;
    varout.ssim_u     = ssim_u;
end
end