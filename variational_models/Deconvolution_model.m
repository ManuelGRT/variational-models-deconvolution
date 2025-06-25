%% --------------------------------------------------------------------- %%
%            Deconvolution: E. Schiavi , I. Ramirez, M. Ramírez           %
%%-----------------------------------------------------------------------%%
function [varout]= Deconvolution_model(varin)
%-------------------------------------------------------------------------%
%                                 Inputs                                  %
%-------------------------------------------------------------------------%
im_org  = varin.im_org;
f       = varin.f;
kernel  = varin.kernel;
Nit     = varin.Nit;
dt      = varin.dt;
epsilon = varin.epsilon_stop;
Verbose = varin.Verbose;

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
kernel_Fourier = psf2otf(kernel,[dim(1),dim(2)]);

%------------------------ Convolution Operations -------------------------%
R  =@(x) real(ifft2(kernel_Fourier.*fft2(x)));
RT =@(x) real(ifft2(conj(kernel_Fourier).*fft2(x)));

%----------------------- Variables initialization ------------------------%
u    = f;
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
            [en(iter),pr(iter),fi(iter)] = pEnergy_R(u,f,1,2,R);
            psnr(iter) = PSNR(u,im_org);
            [ssim_u(iter),ssimmap] = ssim(u,im_org);
            disp(['Iter: ',num2str(iter),' PSNR: ', num2str(psnr(iter)),' SSIM: ', num2str(ssim_u(iter)),' Total Energy: ', num2str(en(iter)), ' Prior: ', num2str(pr(iter)), ' Fidelity: ',num2str(fi(iter)) ])
        case 2
            [en(iter),pr(iter),fi(iter)] = pEnergy_R(u,f,1,2,R);
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
    % Image update
    u2  = u + dt*(-(RT(R(u)-f)));

    % Stop criteria
    norm_u = norm(u2(:) - u(:));
    disp(['Stop Criteria: ',num2str(norm_u)])
    stop = norm_u < epsilon;

    % Variables update
    u = u2;
    iter = iter +1;
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
    varout.ssim_u   = ssim_u;
end
end

