function [u2, L2] = propFF(u1, L1, lambda, z);
% Fraunhofer pattern

[M, N] = size(u1);
dx1 = L1/M;
k = 2*pi/lambda;

L2 = lambda*z/dx1;
dx2 = lambda*z/L1;
x2 = -L2/2 : dx2 : L2/2 - dx2;
[X2, Y2] = meshgrid(x2, x2);

c = 1/(1i*lambda*z) * exp(1i*k/(2*z) * (X2.^2 + Y2.^2));
u2 = c.*ifftshift(fft2(fftshift(u1))) * dx1^2;
end