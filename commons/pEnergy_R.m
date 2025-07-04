function [energy,fidelity,prior]=pEnergy_R(u,f,lambda,p,R)
% Calculates TV(u)+lambda*||u-f||^2
dim = size(u);
Omega = dim(1)*dim(2);
ux       =parcial_x(u);
uy       =parcial_y(u);
modGrad= sqrt(ux.^2+uy.^2 ).^p;
fidelity = (R(u)-f).^2;
prior=sum(modGrad(:))/Omega;
fidelity=sum(fidelity(:))/Omega;
energy= (1/p)*prior+(lambda/2).*fidelity;

end