function [u2] = propIR(u1, L, lambda, z);

[M, N] = size(u1);
dx = L / M;
k = 2*pi / lambda;

x = -L/2 : dx : L/2-dx;
[X, Y] = meshgrid(x, x);

h = 1 / (1i*lambda*z) * exp(1i*k/(2*z)*(X.^2+Y.^2));
H = fft2(fftshift(h)) * dx^2; % dx^2 scaling
U1 = fft2(fftshift(u1));
U2 = H.*U1;
u2 = ifftshift(ifft2(U2));
end