function plap = lap_p(u,p,epsilon)
    % Image Gradient
    ux = parcial_x(u);
    uy = parcial_y(u);

    % Modulus of the Image Gradient
    mod_grad = sqrt(ux.^2 + uy.^2 + epsilon^2);
    b        = mod_grad.^(p-2);

    % Image Laplacian   
    plap = div_x(b.*ux) + div_y(b.*uy);

end