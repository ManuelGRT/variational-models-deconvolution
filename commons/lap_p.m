function plap = lap_p(u,p,epsilon)
    % Image Gradient
    ux = parcial_x(u);
    uy = parcial_y(u);

    % Modulus of the Image Gradient
    mod_grad = sqrt(ux.^2 + uy.^2 + epsilon^2);
    b        = mod_grad.^(p-2);

    % Image Laplacian   
    plap = div_x(b.*ux) + div_y(b.*uy);

    %{
    ux_forward = f([2:end end],:,:)-f;
    uy_forward = f(:,[2:end end],:)-f;
    ux_backward = f-f([1 1:end-1],:,:);
    uy_backward = f-f(:,[1 1:end-1],:);
    ux_mixd = f([2:end end],[1 1:end-1],:)-f(:,[1 1:end-1],:);
    uy_mixd = f([1 1:end-1],[2:end end],:)-f([1 1:end-1],:,:);
    
    plap = zeros(size(f));
    for cc=1:size(f,3)
        plap(:,:,cc)  = ...
            (ux_forward(:,:,cc) + uy_forward(:,:,cc))./max(epsilon,sqrt(ux_forward(:,:,cc).^2 + uy_forward(:,:,cc).^2) ).^(2-p)...
            -ux_backward(:,:,cc)                     ./max(epsilon,sqrt(ux_backward(:,:,cc).^2+ uy_mixd(:,:,cc).^2)    ).^(2-p)...
            -uy_backward(:,:,cc)                     ./max(epsilon,sqrt(ux_mixd(:,:,cc).^2    + uy_backward(:,:,cc).^2)).^(2-p);
    end
    %}
end