function [uout] = tilt(uin, L, lambda, alpha, theta)

[M, N] = size(uin);
dx = L / M;
k = 2*pi/lambda;

x = -L/2 : dx : L/2-dx;
[X, Y] = meshgrid(x, x);

uout = uin.*exp(1i*k*(X*cos(theta) + Y*sin(theta))...
    *tan(alpha));
end