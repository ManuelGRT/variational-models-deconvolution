function ux = parcial_x(u)
% parcial en x forward. es neumann homogenea  ------------------%
ux=zeros(size(u));
ux(:,2:end-1, :)=u(:,3:end,:)-u(:,2:end-1,:);
%-------------------------------------------------------------------------%
end