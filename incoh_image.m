close all;

A = imread('USAF1951w.jpg', 'jpg');
A = imresize(A, [250 250]);
A = rgb2gray(A);

[M, N] = size(A);
A = flipud(A); % reverse row order
Ig = single(A); % integer to floating
Ig = Ig / max(max(Ig)); % normalize ideal image

L = 0.3e-3;
du = L/M;
u = -L/2 : du : L/2-du; v = u;

lambda = 0.5*10^-6;
wxp = 6.25e-3;
zxp = 125e-3;
f0 = wxp / (lambda*zxp);

fu = -1/(2*du) : 1/L : 1/(2*du) - (1/L);
fv = fu;
[Fu, Fv] = meshgrid(fu, fv);
H = circ(sqrt(Fu.^2 + Fv.^2) / f0);

OTF = ifft2(abs(fft2(fftshift(H))).^2); % ?
OTF = OTF / OTF(1, 1);

figure(2)
surf(fu, fv, fftshift(abs(OTF)))
camlight left; lighting phong;
colormap('gray');
shading interp;

Gg = fft2(fftshift(Ig));
Gi = Gg.*OTF;
Ii = ifftshift(ifft2(Gi));
Ii = real(Ii); mask = Ii >= 0; Ii = mask.*Ii;

figure(3)
imagesc(u, v, nthroot(Ii, 2));
colormap('gray');
axis square; axis xy;
