lambda = 0.5e-6;
f = 0.5;
P = 1e-4; % grating period
D1 = 1e-3; % grating side length

L1 = 1e-2;
M = 1000;
dx1 = L1/M;
x1 = -L1/2 : dx1 : L1/2 -dx1;
[X1, Y1] = meshgrid(x1, x1);
