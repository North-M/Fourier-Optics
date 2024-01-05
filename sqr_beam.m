close all;

L1 = 0.5;
M = 250;
dx = L1/M;
dx1 = L1/M;
x1 = -L1/2 : dx1 : L1/2-dx1;
y1 = x1;

lambda = 0.5*10^-6;
k = 2*pi/lambda;
w = 0.051;
z = 20000;

[X1, Y1] = meshgrid(x1, y1);
u1 = rect(X1 / (2*w)).*rect(Y1/(2*w));
I1 = abs(u1.^2);

figure(1);
imagesc(x1, y1, I1);
axis square; axis xy;
colormap('gray');

u2 = propTF(u1, L1, lambda, z);
% u2 = propIR(u1, L1, lambda, z);

x2 = x1; y2 = y1;
I2 = abs(u2.^2);

figure(2);
imagesc(x2, y2, I2);
axis square; axis xy;
colormap('gray');

figure(3);
plot(x2, I2(M/2+1, :));

figure(4);
plot(x2, abs(u2(M/2+1, :)));

figure(5);
plot(x2, unwrap(angle(u2(M/2+1, :))));

w = 0.011;
u1 = rect(X1 / (2*w)).*rect(Y1/(2*w));
z = 200000;
[u2, L2] = propFF(u1, L1, lambda, z);

dx2 = L2/M;
x2 = -L2/2 : dx2 : L2/2-dx2;
y2 = x2;
I2 = abs(u2.^2);

figure(6);
imagesc(x2, y2, nthroot(I2, 3))
axis square; axis xy;
colormap('gray');
