% I = rgb2gray(imread('sample.jpg'));
% Sigma=estimate_noise(I);
%         
% The advantage of this method is that it includes a Laplacian operation 
% which is almost insensitive to image structure but only depends on the 
% noise in the image. 

function Sigma=estimate_noise_rgb(I)

    [H,W,C]=size(I);
    I=double(I);
    temp = [];
    % compute sum of absolute values of Laplacian
    M=[1 -2 1; -2 4 -2; 1 -2 1];
    %M=[-1 -1 -1; -1 8 -1; -1 -1 -1];
    %M=[0 -1 0; -1 4 -1; 0 -1 0];

    if C > 1
        for i=1:C
            temp(i)=sum(sum(abs(conv2(I(:,:,i), M))));
        end
        Sigma = sum(temp)/C;
    else
        Sigma=sum(sum(abs(conv2(I, M))));
    end

    % scale sigma with proposed coefficients
    Sigma=Sigma*sqrt(0.5*pi)./(6*(W-2)*(H-2));

end