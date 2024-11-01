function d = loadDoe(diffractive,N,sigma_d)
if (diffractive == 0)
    load('Heightmap.mat');
    d = heightmap;
    d = imresize(d,[N N],"nearest");
elseif(diffractive == 1)
    %load('filter1800_16.mat');
    %load('filter_big_31800_16.mat')
    %load('filter_big_61800_16.mat')
    %load('filter_big_71800_16.mat')
    %load('filter1800_16_3.mat');
    %load('filter1800_16_2.mat')
    %load('filter1800_16_few.mat')
    %load('filter1800_16_one.mat')
    %load('Uniform-DOE_1800.mat')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spiral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load('Uniform-DOE_1800-clas-0.5.mat')
    M1 = 2464;
    load('Uniform-DOE_1800-fres-0.5-0.8.mat')
    D1 = imresize(G3,[M1 M1],"nearest");
    load('Uniform-DOE_1800-clas-0.8-1.2.mat')
    D2 = imresize(G3,[M1 M1],"nearest");
    load('Uniform-DOE_1800_0.8-1.2-inf.mat')
    D3 = imresize(G3,[M1 M1],"nearest");
    a = 0.1;
    b = 0.09;
    c = 1-a-b;
    
    im = ones(M1,M1);
    A_st1 = elliptical_crop(im,0.5)>0;
    A_st2 = elliptical_crop(im,0.9)>0;
    A_st3 = elliptical_crop(im,1)>0;
    %A_st4 = elliptical_crop(im,0.8)>0;
    % a = 0.45;
    % b = 0.05;
    % c = 0.65;
    % + b*D2 + c*D3; 
    d =     1*D1.*(A_st3-A_st2) + 1*D1.*A_st1 + 1*(A_st2-A_st1).*D3 ;
    %A_st4 = elliptical_crop(im,(3/9.2))>0;
    %d = d.*A_st4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pyramid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load('Uniform-DOE_1800_0.5.mat')
    % D1 = G3;
    % load('Uniform-DOE_1800_0.8-1.2-inf.mat')
    % D2 = G3;
    % a = 0.4;
    % b = 0.6;
    % G3 = a*D1 + b*D2;


    %9.2 ->2464
    %3 -> x
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    d = imresize(d,[N N],"nearest");
    d = d.*1.3745e-06;%.*1.1000e-06; %.*1.3745e-06;%.*1e-06; % 1.1348e-06;
    G3 = d;
    save('./DOEs/Uniform-DOE_1800.mat','G3')
elseif(diffractive == 2)
    %Paper: Compact Snapshot Hyperspectral Imaging with Diffracted Rotation
    %Authors: Jeon, Baek, Yi, Dun, Fu, Heidrich, and Kim
    load('spiral_DOE3180_16.mat');
    d = h;
    d = imresize(d,[N N],"nearest");
    d = mat2gray(d).*0.9e-06;
elseif(diffractive == 3)
    alpha = 10*pi; % 5 10 -> inf % 10 d= 1-d -> 0.5 0.8
    M = 1800;
    delta = 3e-6; %pixel size
    [x, y] = meshgrid((-M/2 : M/2-1) *pi*delta);
    d = 1*(x./max(x(:))).^2 + 1*(y./max(y(:))).^2;
    d = mat2gray(mod(alpha.*d, 2*pi));
    d = d;
    %d = d.*1.3745e-06;
    im = ones(M,M);
    A_st = elliptical_crop(im,1)>0;
    d = (1-d);
    imagesc(d)
    %load('fresnel3180_16.mat');
    %d = G3;
    d = imresize(d,[N N],"nearest");
    d = mat2gray(d).*1.3745e-06;%.*0.9e-06;
    %G3 = d;
    save('./DOEs/fresnel-DOE_1800.mat','G3')
    %load('Uniform-DOE_1800-fres-0.5-0.8.mat')
elseif(diffractive == 4)
    d = zeros(N,N);
elseif(diffractive == 5)
    %Paper: Extended depth of field through wave-front coding
    % %AuthorS: Edward R. Dowski, Jr., and W. Thomas Cathey
    alpha = 20*pi;
    M = N;
    delta = 3e-6; %pixel size
    [x, y] = meshgrid((-N/2 : M/2-1) *pi*delta);
    d = 1*(x./max(x(:))).^3 + 1*(y./max(y(:))).^3;
    d = mat2gray(mod(alpha.*d, 2*pi));
    d = d.*1.3745e-06;
    im = ones(N,N);
    A_st = elliptical_crop(im,1)>0;
    d = d.*A_st;
    save('./DOEs/Dowski-Cathey_test.mat','d')
    %load('Dowski-Cathey_test.mat')
elseif(diffractive==6)
    %load('filter1800_16.mat');
    %load('filter_big_31800_16.mat')
    %load('filter_big_61800_16.mat')
    %load('filter_big_71800_16.mat')
    %load('filter1800_16_3.mat');
    %load('filter1800_16_2.mat')
    %load('filter1800_16_few.mat')
    %load('filter1800_16_one.mat')
    %load('Uniform-DOE_1800.mat')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spiral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load('Uniform-DOE_1800-clas-0.5.mat')
    load('Uniform-DOE_1800-fres-0.5-0.8.mat')
    D1 = G3;
    %load('Uniform-DOE_1800-clas-0.8-1.2.mat')
    load('Uniform-DOE_1800-clas-0.8.mat')
    D2 = G3;
    load('Uniform-DOE_1800_0.8-1.2-inf.mat')
    D3 = G3;
    % a = 0.1;
    % b = 0.09;
    % c = 1-a-b;
    M1 = 1800;
    im = ones(M1,M1);
    A_st1 = elliptical_crop(im,0.2)>0;
    A_st2 = elliptical_crop(im,0.65)>0;
    A_st3 = elliptical_crop(im,1)>0;
    %A_st4 = elliptical_crop(im,0.8)>0;
    % a = 0.45;
    % b = 0.05;
    % c = 0.65;
    % + b*D2 + c*D3; 
    d =     1*(D1).*(A_st3-A_st2) + (1*(D2).*A_st1) + 1*(A_st2-A_st1).*(D3) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pyramid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load('Uniform-DOE_1800_0.5.mat')
    % D1 = G3;
    % load('Uniform-DOE_1800_0.8-1.2-inf.mat')
    % D2 = G3;
    % a = 0.4;
    % b = 0.6;
    % G3 = a*D1 + b*D2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    d = imresize(d,[N N],"nearest");
    d = d.*1.3745e-06;%.*1.1000e-06; %.*1.3745e-06;%.*1e-06; % 1.1348e-06;
    G3 = d;
    save('./DOEs/Uniform-DOE_1800.mat','G3')
elseif(diffractive==7)    
    load('Phase_delay_profile_6mm.mat')  % no full
    G3 = Phase_delay_profile_6mm_crop_middle;
    % load('Wrapped_Phase_delay_profile_6mm_0_2pi.mat'); % no full
    % G3 = Wrapped_Phase_delay_profile_6mm;
    load('Phase_delay_profile_9.2mm.mat') % full
    G3 = Phase_delay_profile;
    % load('Wrapped_Phase_delay_profile_9.2mm_0_2pi.mat'); % full
    % G3 = Wrapped_Phase_delay_profile_9mm; 
    length(unique(G3(:)))
    max(G3(:))
    d = imresize(G3,[N N],"nearest");
    d = mat2gray(d).*1.3745e-06;
    %imagesc(d)
    
end
if(diffractive~=4)
    d = d + normrnd(0,sigma_d,[N,N]);
end