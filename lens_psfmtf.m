close all;

M = 1024;
L = 1e-3;
du = L / M;
u = -L/2 : du : L/2 - du; v = u;

lambda = 0.55e-6;
k = 2*pi / lambda;
Dxp = 20e-3; wxp = Dxp / 2;
zxp = 100e-3;

fnum = zxp / (2*wxp);
lz = lambda*zxp;
twof0 = 1 / (lambda*fnum);

u0 = 0; v0 = 0.5;

wd = 0;
w040 = 4.963*lambda;
w131 = 2.637*lambda;
w222 = 9.025*lambda;
w220 = 7.536*lambda;
w311 - 0.157*lambda;

fu = -1/(2*du) : 1/L : 1/(2*du) - (1/L);
[Fu, Fv] = meshgrid(fu, fu);

% wavefront
W = seidel_5(u0, v0, -lz*Fu/wxp, -lz*Fv/wxp, wd, w040, w131, w222, w220, w311);

% coherent transfer function
H = circ(sqrt(Fu.^2 + Fv.^2)*lz/wxp).*exp(-1i*k*W);

figure(1)
imagesc(u, v, angle(H));
axis xy; axis square;
colormap('gray');

% point spread function
h2 = abs(ifftshift(ifft2(fftshift(H)))).^2;

figure(2)
imagesc(u, v, nthroot(h2, 2));
axis xy; axis square;
colormap('gray');

figure(3)
plot(u, h2(M/2+1, :));

figure(4)
plot(u, h2(:, M/2+1));

% MTF
MTF = fft2(fftshift(h2));
MTF = abs(MTF/MTF(1, 1));
MTF = ifftshift(MTF);

MTF_an = (2/pi) * (acos(fu/twof0) - (fu/twof0).*sqrt(1 - (fu/twof0).^2));
MTF_an = MTF_an.*rect(fu/(2*twof0));

figure(5)
plot(fu, MTF(M/2+1, :), fu, MTF(:, M/2+1), ':', fu, MTF_an, '--');
axis([0 150000 0 1]);
legend('u MTF', 'v MTF', 'diff limit');
