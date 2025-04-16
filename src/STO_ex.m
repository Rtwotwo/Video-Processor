clc;
close all;
clear all;

Nc = 64;

S_r = 0.5^0.5 * ((randn(1, Nc) > 0) * 2 - 1);
S_i = 0.5^0.5 * ((randn(1, Nc) > 0) * 2 - 1);

S = S_r + 1i * S_i;

STO = -23;

y = S .* exp(1i * 2 * pi * STO / Nc * (0 : Nc-1));

figure()
plot(S_r, S_i, '.');

figure()
plot(real(y), imag(y), '.');