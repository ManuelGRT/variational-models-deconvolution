%% --------------------------------------------------------------------- %%
%     Deconvolution and Denoising: E. Schiavi, I. Ramírez, M. Ramirez     %
%%-----------------------------------------------------------------------%%
function [varout]= Deconvolution_Denoising_model(varin)
%-------------------------------------------------------------------------%
%                                 Inputs                                  %
%-------------------------------------------------------------------------%
im_org   = varin.im_org;
f        = varin.f;
kernel   = varin.kernel;

p        = varin.p;
Nit      = varin.Nit;
dt       = varin.dt;
lambda1  = varin.lambda1;
lambda2  = varin.lambda2;
epsilon  = varin.epsilon;
hasNoise = varin.hasNoise;
stop_C   = varin.stop_C;
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

%---------------- Kernel Transformation to Fourier Space -----------------%
dim      = size(f);
kernel_F = psf2otf(kernel,[dim(1),dim(2)]);

%------------------------ Convolution Operations -------------------------%
R  =@(x) real(ifft2(kernel_F.*fft2(x)));
RT =@(x) real(ifft2(conj(kernel_F).*fft2(x)));

%----------------------- Variables initialization ------------------------%
u    = f;
iter = 1;
stop = false;
if hasNoise
    sigma_hat = estimate_noise_rgb(f);
end

while (iter <= Nit) && (~stop)
%---------------------------- Visualization ------------------------------%
    switch Verbose
        case 0

            disp(['Iter: ',num2str(iter)])
            
        case 1
            if ~mod(iter,100)
                imshow(u),pause(0.01)
            end
            [en(iter),pr(iter),fi(iter)] = pEnergy_R(u,f,lambda1,p,R);
            psnr(iter) = PSNR(u,im_org);
            [ssim_u(iter),ssimmap] = ssim(u,im_org);
            disp(['Iter: ',num2str(iter),' PSNR: ', num2str(psnr(iter)),' SSIM: ', num2str(ssim_u(iter)),' Total Energy: ', num2str(en(iter)), ' Prior: ', num2str(pr(iter)), ' Fidelity: ',num2str(fi(iter)) ])
        case 2
            [en(iter),pr(iter),fi(iter)] = pEnergy_R(u,f,lambda1,p,R);
            psnr(iter) = PSNR(u,im_org);
            [ssim_u(iter),ssimmap] = ssim(u,im_org);
            subplot(141), imshow(u), title(['PSNR: ',num2str(psnr(iter))])
            subplot(142), plot(en,'r','LineWidth', 2), hold on,
                          plot(fi,'b','LineWidth', 2),
                          plot(pr,'g','LineWidth', 2)
                          legend('Total Energy','Fidelity','Prior'), grid on
            subplot(143), plot(psnr,'c','LineWidth', 2), 
                          legend('PSNR (db)'), grid on
            subplot(144), imshow(ssimmap,[]), title("SSIM: "+ssim_u(iter))
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

    % Image update by Gradient Descense
    mask = exp(mod_grad);

    u2  = u - dt*(mask.*(RT(R(u)-f)) - plap);
    %u2  = u - dt*(lamda1*(RT(R(u)-f)) - lambda2*plap);

    if hasNoise 
        % Stop criteria 1: Image with noise
        size_omega = numel(f);           % Number of pixels
        residual = (f - R(u)).^2;
        integral = sum(residual(:));

        disp(['Stop Criteria: ', num2str(integral - size_omega*stop_C*(sigma_hat.^2))])
        stop = integral - size_omega*stop_C*(sigma_hat.^2) < 0;

    else
        % Stop criteria 2: Image without noise
        norm_u = norm(u2(:) - u(:));
        disp(['Stop Criteria: ', num2str(norm_u)])

        stop = norm_u < epsilon;    
    end

    % Variables update
    u = u2;
    iter = iter +1;

    %lambda1 = max(10,1.01*lambda1);
    %lambda2 = max(1e-3,0.99*lambda2);
end

%--------------------------- Output Variables ----------------------------%
if Verbose == 0
    varout.u        = u;
else
    varout.u        = u;
    varout.en       = en;
    varout.pr       = pr;
    varout.fi       = fi;
    varout.psnr     = psnr;
    varout.ssim_u     = ssim_u;
end
end

