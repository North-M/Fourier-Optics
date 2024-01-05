close all;

A = imread('USAF1951w.jpg', 'jpg');
A = imresize(A, [250 250]);
A = rgb2gray(A);

[M, N] = size(A);
A = flipud(A); % reverse row order
Ig = single(A); % integer to floating
Ig = Ig / max(max(Ig)); % normalize ideal image

% ug = sqrt(Ig); % ideal image field;
ug = sqrt(Ig).*exp(1i*2*pi*rand(M)/2);
L = 0.3e-3;
du = L/M;
u = -L/2 : du : L/2-du; v = u;

figure(1)
imagesc(u, v, Ig);
colormap('gray');
axis square; axis xy;

lambda = 0.5*10^-6;
wxp = 1.25e-2;
zxp = 125e-3;
f0 = wxp / (lambda*zxp);

fu = -1/(2*du) : 1/L : 1/(2*du) - (1/L);
fv = fu;
[Fu, Fv] = meshgrid(fu, fv);
H = circ(sqrt(Fu.^2 + Fv.^2) / f0);

figure(2)
surf(fu, fv, H.*.99);
camlight left; lighting phong;
colormap('gray');
shading interp;

H = fftshift(H);
Gg = fft2(fftshift(ug));
Gi = Gg.*H;
ui = ifftshift(ifft2(Gi));
Ii = (abs(ui)).^2;

figure(3)
imagesc(u, v, nthroot(Ii, 2));
colormap('gray');
axis square; axis xy;

figure(4)
vvalue = -0.8e-4;
vindex = round(vvalue/du + (M/2 + 1));
plot(u, Ii(vindex, :), u, Ig(vindex, :), ':');

