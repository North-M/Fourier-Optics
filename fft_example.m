% fft_example

close all;

w = 0.055;
L = 2;
M = 200;
dx = L / M;

x = -L/2:dx:L/2-dx;
f = rect(x / (2*w));

figure(1);
plot(x, f);

figure(2);
plot(x, f, '-o');
axis([-0.2 0.2 0 1.5]);
xlabel('x(m)');

figure(3);
plot(f, '-o');
axis([80 120 0 1.5]);
xlabel('index');

f0 = fftshift(f);
figure(4);
plot(f0);
axis([0 200 0 1.5]);
xlabel('index');

F0 = fft(f0) * dx;
figure(5);
plot(abs(F0));
title('magnitude');
xlabel('index');

figure(6);
plot(angle(F0));
title('phase');
xlabel('index');

F = fftshift(F0);
fx = -1/(2*dx):1/L:1/(2*dx)- 1/L;
figure(7);
plot(fx, abs(F));
title('magnitude');
xlabel('fx cyc/m');

figure(8);
plot(fx, angle(F));
title('phase');
xlabel('fx cyc/m');

F_an = 2*w*sinc(2*w*fx);

figure(9);
plot(fx, abs(F), fx, abs(F_an), ':');
title('magnitude');
legend('discrete', 'analytic');

figure(10);
plot(fx, angle(F), fx, angle(F_an), ':');
title('phase');
legend('discrete', 'analytic');
