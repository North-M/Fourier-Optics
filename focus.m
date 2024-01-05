function [uout] = focus(uin, L, lambda, zf)

[M, N] = size(uin);
dx = L/M;
k = 2*pi/lambda;

x = -L/2 : dx : L/2-dx;
[X, Y] = meshgrid(x, x);

uout = uin.*exp(-1i*k/(2*zf) * (X.^2 + Y.^2));
end