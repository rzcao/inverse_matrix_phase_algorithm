function [ noiseimage ] = quasi_wnr( OTF,noiseimage,noise_factor)
%low-pass filter
[~,~,num]=size(noiseimage);
abs_OTF=abs(OTF);
conj_OTF=conj(OTF);
for ii=1:num
    ft=fftshift(fft2(noiseimage(:,:,ii)));
    filtered_ft=(ft.*conj_OTF./(noise_factor+abs_OTF));
    noiseimage(:,:,ii)=ifft2(ifftshift(filtered_ft));
end

end

