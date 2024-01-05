L = 0.2; % side length (m)
M = 250; % # samples
dx = L/M; % sample interval
x = -L/2 : dx : L/2-dx; y = x; % corrds
[X, Y] = meshgrid(x, y);

w = 1e-3; % x half-width of aperture
lambda=0.633e-6; % wavelength (um)
z = 50; % prop distance (m)
k = 2*pi/lambda; % wavenumber
lz = lambda*z;

%irradiance
I2 = (w^2 / lz)^2.*(jinc(w/lz * sqrt(X.^2 + Y.^2))).^2;

figure(1)
imagesc(x, y, nthroot(I2, 3));
% imagesc(x, y, I2);
colormap('gray');
axis square;
axis xy;

figure(2)
plot(x, I2(M/2+1, :));