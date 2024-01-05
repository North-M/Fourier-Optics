close all;

M = 250;
L = 1e-3;
du = L / M;
u = -L/2 : du : L/2-du; v = u;

lambda = 0.5e-6;
k = 2*pi/lambda;
wxp = 2.5e-3;
zxp = 100e-3;
fnum = zxp / (2*wxp);

twof0 = 1 / (lambda*fnum);
fN = 1 / (2*du);

wd = 0;
w040 = 0.5*lambda;
w131 = 1*lambda;
w222 = 1.5*lambda;
w220 = 0;
w311 = 0;

fu = -1/(2*du) : 1/L : 1/(2*du) - (1/L);
fu = fftshift(fu); % shift cords, avoid shifting H in loop
[Fu, Fv] = meshgrid(fu, fu);

I = zeros(M);
for u0 = [-0.7:0.7/3:0.7]
    for v0 = [-0.7:0.7/3:0.7]
        % wavefornt
        W = seidel_5(u0, v0, -2*lambda*fnum*Fu...
            , -2*lambda*fnum*Fv, ...
            wd, w040, w131, w222, w220, w311);
        % coherent transfer function
        H = circ(sqrt(Fu.^2 + Fv.^2) * lambda * fnum)...
            .*exp(-1i*k*W);
        % PSF
        h2 = abs(ifftshift(ifft2(H))).^2;
        % shift PSF to image plane position
        h2 = circshift(h2, [round(v0*M/2), round(u0*M/2)]);
        I = I + h2;
    end
end

figure(1)
imagesc(u, v, nthroot(I, 2));
colormap('gray'); axis square; axis xy;