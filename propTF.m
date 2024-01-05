function [u2] = propTF(u1, L, lambda, z);

[M, N] = size(u1);
dx = L / M;
k = 2*pi / lambda;

fx = -1/(2*dx) : 1/L : 1/(2*dx) - 1/L;
[FX, FY] = meshgrid(fx, fx);

H = exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));
H = fftshift(H);
U1 = fft2(fftshift(u1));
U2 = H.*U1;
u2 = ifftshift(ifft2(U2));
end