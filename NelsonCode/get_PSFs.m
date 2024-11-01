function [rgbPSF,rgbPSF_mis,psf_infocus,psf_misfocus] = get_PSFs(x0,phi,specH)

LevelNum  = x0(1);
FrenOrder = x0(2);

%%

lambda_0   = 5.1e-7;
image_size = 512;
ny_sensor  = image_size; nx_sensor=image_size;
ny = ny_sensor/1; nx = nx_sensor/1;
pixel_size = 3.45e-6;
delta_u    = pixel_size;

f_0 = 0.01;
D   = 0.0092;
a22 = 1e-2;
z1_focus = 1;
n_0 = OpticalMedium('BK7').refractiveIndex(lambda_0*10^9, 'nm');

z1var = [0.5 1.0 1.2];
z2    = (1/f_0 - 1/z1_focus)^(-1);


N       = 2464;

delta_x = D/N;

ux = linspace(-delta_u*nx/2, delta_u*nx/2,nx);
uy = linspace(-delta_u*ny/2, delta_u*ny/2,ny);
x  = linspace(-delta_x*N/2,delta_x*N/2,N);

y=x;

xrow = full(x(:)).';
ycol = full(y(:));
X = repmat(xrow,size(ycol));
Y = repmat(ycol,size(xrow));

wavelengths = 1e-9*(400:10:700);

weights = [0.1011 0.08079 0.0587 0.04371 0.0366 0.03046 0.02517 0.02693 0.03488 0.04283 0.04812 0.05695 0.07285 0.08874 0.08256 0.06402 0.0578 0.0755 0.1912 0.4667 0.5779 0.5532 0.5311 0.5196 0.4923 0.4517 0.4102 0.3404 0.2141 0.09669 0.04194;
    0.102 0.0896 0.0763 0.064 0.06578 0.07285 0.07815 0.1488 0.3333 0.5929 0.7166 0.7563 0.7607 0.7563 0.7298 0.6848 0.6327 0.5744 0.5143 0.4437 0.3581 0.2565 0.1779 0.1329 0.1073 0.09404 0.0843 0.07991 0.06049 0.03488 0.02075;
    0.3342 0.45 0.5205 0.5859 0.6362 0.6768 0.683 0.6583 0.6283 0.5558 0.4517 0.3395 0.2389 0.1629 0.1249 0.0931 0.0657 0.0481 0.0437 0.04106 0.0357 0.0278 0.02428 0.02517 0.02693 0.03046 0.03488 0.03664 0.02781 0.01545 0.0083];

for ii=1:3
    norm_weights(ii,:) = weights(ii,:)./sum(weights(ii,:));
end


R     = 3e-3;
pupil = (X.^2 + Y.^2 <= R.^2);
size_wavelength = size(wavelengths);

[~, size_z1] = size (z1var);


psf_infocus  = zeros(size_wavelength(2),image_size,image_size);
psf_misfocus = zeros(size_z1,size_wavelength(2),image_size,image_size);

%%
for i=1:size_wavelength(2)


    lambda  = wavelengths(i);
    n       = OpticalMedium('BK7').refractiveIndex(lambda*10^9, 'nm');

    specH = (2*pi)*(n_0-1)/lambda_0*specH;
    phi     = zeros(size(phi));

    f = f_0*(n_0-1)/(n-1);

    phi_wrapped_Fresnel_Order = mod(phi+pi*FrenOrder,2*(pi*FrenOrder))-(pi*FrenOrder);
    mask = round( phi_wrapped_Fresnel_Order/(2*pi*FrenOrder/LevelNum) )*(2*pi*FrenOrder/LevelNum);

    BPM_0 = pupil.*mask;

    h     = BPM_0*lambda_0/(n_0-1)/2/pi;

    BPM_corrected = (2*pi)*(n-1)/lambda*h;
    %BPM_corrected = (BPM_corrected + abs(min(BPM_corrected(:)))).*pupil;

    for delta = 1:size_z1

        z1 = z1var(delta);

        if z1 == z1_focus

            Pg_0  = pupil.*exp(1i*pi/lambda*(1/z1+1/z2-(1-a22)/f).*(X.^2 + Y.^2) + 1i*BPM_corrected);
            XX    = diag(sin(pi/lambda/z2*x*delta_u)./(pi/lambda/z2*x));
            A     = exp(-1j*2*pi*uy'*x/lambda/z2)*XX; B=XX*exp(-1j*2*pi*x'*ux/lambda/z2);
            C_PSF = A*Pg_0*B;

            PSF   = abs(C_PSF).^2;
            psf_infocus(i,:,:) = PSF/sum(PSF(:));
        end

        Pg_0 = pupil.*exp(1i*pi/lambda*(1/z1+1/z2-(1-a22)/f).*(X.^2 + Y.^2) + 1i*BPM_corrected);

        XX    = diag(sin(pi/lambda/z2*x*delta_u)./(pi/lambda/z2*x));
        A     = exp(-1j*2*pi*uy'*x/lambda/z2)*XX; B=XX*exp(-1j*2*pi*x'*ux/lambda/z2);
        C_PSF = A*Pg_0*B;

        PSF   = abs(C_PSF).^2;

        psf_misfocus(delta,i,:,:) = PSF/sum(PSF(:));

    end
end

rgbPSF = zeros(3,image_size,image_size);

for i=1:3
    for jj=1:size_wavelength(2)
        rgbPSF(i,:,:)=rgbPSF(i,:,:) + psf_infocus(jj,:,:)*norm_weights(i,jj);
    end
    tmp = rgbPSF(i,:,:);
    rgbPSF(i,:,:) = rgbPSF(i,:,:)/sum(tmp(:));
end

rgbPSF_mis = zeros(size_z1,3,image_size,image_size);

for delta = 1:size_z1
    for i=1:3
        for jj=1:31
            rgbPSF_mis(delta,i,:,:) = rgbPSF_mis(delta,i,:,:) + psf_misfocus(delta,jj,:,:)*norm_weights(i,jj);
        end
        tmp = rgbPSF_mis(delta,i,:,:);
        rgbPSF_mis(delta,i,:,:) = rgbPSF_mis(delta,i,:,:)/sum(tmp(:));
    end
end
end


