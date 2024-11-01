clear;
close all
clc


%% optimized

levelnum = 52;
frenorder = 1.4;


x0 = [levelnum, frenorder];
load('Phase_delay_profile_6mm.mat','Phase_delay_profile_6mm_crop_middle');
load('Uniform-DOE_2464-.mat','G3');

phi = Phase_delay_profile_6mm_crop_middle;
specH = G3;

[rgbPSF,rgbPSF_mis,psf_wav,psf_wav_miss] = get_PSFs(x0,phi,specH);

%%

wavelengths = 1e-9*(400:10:700);

figure;
for l=1:size(psf_wav,1)
    imagesc(squeeze(psf_wav(l,:,:))),title(num2str(wavelengths(l)))
    pause(0.2)
end
    
