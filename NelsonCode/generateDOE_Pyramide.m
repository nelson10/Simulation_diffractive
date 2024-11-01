clear;
clc;
close all;
%% Add path to the dataset

N  = 1600;
NF = N/2;

NF1 = NF;
A = zeros(NF1,NF1);
for i=1:NF1
    A(i,:) = circshift(NF1:-1:1,-(i-1)); % Uniform (Aggarwal)
end

A = imrotate(A,0);
A2 = imrotate(A,-270);
A3 = imrotate(A,-180);
A1 = imrotate(A,270);
A = [A3 A2;A1 A];
A = imresize(A,[NF,NF]);

p = ceil(N/NF1);
B1 = ones(p,p);

imagesc(A)
G = kron(B1,A);
G = G(1:N,1:N);

img = G/max(G(:));
[h,w,~] = size(img);
s = min(h,w)/2;
[rho,theta] = meshgrid(linspace(0,s-1,s), linspace(0,2*pi));
[x,y] = pol2cart(theta, rho);
z = zeros(size(x));
figure(1)
subplot(121), imshow(img)
subplot(122),warp(x, y, z, img), view(2), axis square tight off
F = getframe ;

I = F.cdata;

%%
I = rgb2gray(I);
im = ones(N,N);
ap = elliptical_crop(im,1)>0;
Ir = imresize(I,[N N],method="nearest");
Ir = mat2gray(Ir);
Ir = ap.*double(Ir);

G = Ir;
G3 = G;
G3 = imrotate(G3,-90);

figure;imagesc(G3)

%%

M = 2464;
finalDOE = zeros(M,M);
finalDOE(M/2-N/2:M/2+N/2-1,M/2-N/2:M/2+N/2-1) = G3;

figure;imagesc(finalDOE)
