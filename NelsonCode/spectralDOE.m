clc;
close all;
clear;
%%

alpha = 10*pi;
M     = 1600;
D     = 6e-3;
R     = 3e-3;
delta = D/M;

x     = linspace(-delta*M/2,delta*M/2,M);

[X,Y] = meshgrid(x,x);
d     = 1*(X./max(X(:))).^2 + 1*(Y./max(Y(:))).^2;
d     = mat2gray(mod(alpha.*d, 2*pi))*1.3745e-06;

pupil = double((X.^2 + Y.^2)<=R.^2);
d = d.*pupil;

figure;imagesc(d)

%%
N = 2464;
finalDOE = zeros(N,N);
finalDOE(N/2-M/2:N/2+M/2-1,N/2-M/2:N/2+M/2-1) = d;

figure;imagesc(finalDOE)
