function uy = parcial_y(u)
% parcial en Y forward. Es neumann homogenea --- ---------------%
uy=zeros(size(u));
uy(2:end-1, :,:)=u(3:end,:,:)-u(2:end-1,:,:);
%-------------------------------------------------------------------------%
end