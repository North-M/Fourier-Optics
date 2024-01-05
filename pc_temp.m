close all;

lambda0 = 650e-9; % center wavelength (m)
c = 3e8;
k0 = 2*pi / lambda0;
nu0 = c / lambda0; % center frequency (Hz)

% Gaussian linshape parameters
N = 51;
delnu = 2e9; % spectral desity FWHM (Hz)
b = delnu / (2*sqrt(log(2))); % FWHM scaling
dnu = 4 * delnu / N; % freq interval (?)

% source plane parameters
L1 = 50e-3;
M = 250;
dx1 = L1 / M;
x1 = -L1/2 : dx1 : L1/2-dx1;
x1 = fftshift(x1);
[X1, Y1] = meshgrid(x1, x1);

% beam parameters
w = 1e-3;
dels = 5e-3; % transverse sparation
deld = 5e-2; % delay distance
f = 0.25; % focal dist for Fraunhofer
lf = lambda0 * f;

% loop through lines
I2 = zeros(M);
for n = 1:N
    % soectral density function
    nu = (n - (N+1)/2) * dnu + nu0;
    S = 1/(sqrt(pi)*b) * exp(-(nu-nu0)^2 / b^2);
    k = 2*pi*nu/c;
    % source
    u = circ(sqrt((X1 - dels/2).^2 + Y1.^2) / w)...
        + circ(sqrt((X1+dels/2).^2 + Y1.^2) / w) * exp(1i*k*deld);
    % Fraunhofer pattern
    u2 = 1/lf * (fft2(u)) * dx1^2;
    % weighted irradiance and sum
    I2 = I2 + S * (abs(u2).^2) * dnu;
end

I2 = ifftshift(I2);
x2 = (-1/(2*dx1) : 1/L1 : 1/(2*dx1)-1/L1) * lf; % obs corrds
y2 = x2;

figure(1)
imagesc(x2, y2, nthroot(I2, 3));
axis square; axis xy; colormap('gray');

figure(2)
plot(x2, I2(M/2+1, :));
