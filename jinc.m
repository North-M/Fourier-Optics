function [out] = jinc(x);

% J1(2*pi*x) / x

% locate non-zero elements of x
mask = (x ~= 0);

% initialize output with pi (value for x = 0)

out = pi*ones(size(x));

% compute output values for all other x

out(mask) = besselj(1, 2*pi*x(mask))./x(mask);
end