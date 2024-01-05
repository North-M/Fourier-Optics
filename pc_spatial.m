close all;

lambda = 650e-9;

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
f = 0.25; % focal dist for Fraunhofer
lf = lambda0 * f;

% partial spatial coherence screen parameters
N = 100; % number of screens
Lcr = 8e-3; % spatial correlation length
sigma_f = 2.5*Lcr; % Gaussian filter parameter
sigma_r = sqrt(4*pi*sigma_f^4 / Lcr^2); % std

dfx1 = 1/L1;
fx1 = -1/(2*dx1) : dfx1 : 1/(2*dx1)-dfx1;
fx1 = fftshift(fx1);
[FX1, FY1] = meshgrid(fx1, fx1);

% source field
u1 = circ(sqrt((X1 - dels/2).^2 + Y1.^2) / w)...
        + circ(sqrt((X1+dels/2).^2 + Y1.^2) / w);
% filter spectrum
F = exp(-pi^2*sigma_f^2 * (FX1.^2 + FY1.^2));

% loop through screens
I2 = zeros(M);
for n = 1:N/2
    % make 2 random screens (real and imag)
    fie = (ifft2(F.*(randn(M) + 1i*randn(M)))...
        * sigma_r/dfx1) * M^2 * dfx1^2;
    % Fraunhofer pattern applying screen 1
    u2 = 1/lf * (fft2(u1.*exp(1i*real(fie)))) * dx1^2;
    I2 = I2 + abs(u2).^2;
    % Fraunhofer pattern applying screen 2
    u2 = 1/lf * (fft2(u1.*exp(1i*imag(fie)))) * dx1^2;
    I2 = I2 + abs(u2).^2;
end

I2 = ifftshift(I2) / N;
x2 = (-1/(2*dx1) : 1/L1 : 1/(2*dx1)-1/L1)*lf;
y2 = x2;

figure(1)
imagesc(x2, y2, nthroot(I2, 3));
axis square; axis xy;
colormap('gray');

figure(2)
plot(x2, I2(M/2 + 1, :));
