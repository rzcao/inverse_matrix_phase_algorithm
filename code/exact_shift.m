function [shift_im] = exact_shift(im,relative_pixel,size_flag)
% shift the matrix (subpixel supported)
% similar to the built-in function: circshift
%size_flage=1: output_matrix is the same size as the input matrix
[xsize,ysize,~]=size(im);
if size_flag==1
    o_xsize=xsize;
    o_ysize=ysize;
end

if mod(xsize,2)==0
    im(xsize+1,:,:)=0;
end

if mod(ysize,2)==0
    im(:,ysize+1,:)=0;
end
relative_pixel(:,1)=rem(relative_pixel(:,1),xsize);
relative_pixel(:,2)=rem(relative_pixel(:,2),ysize);
[xsize,ysize,num]=size(im);
[Y,X]=meshgrid(1:ysize,1:xsize);
xc=floor(xsize/2+1);
yc=floor(ysize/2+1);
yr=Y-yc;
xr=X-xc;
r_shift=sqrt(relative_pixel(:,1).^2+relative_pixel(:,2).^2);
shift_angle=zeros(num);
f_r=zeros(num);
for ii=1:num
    if relative_pixel(ii,2)==0
        shift_angle(ii)=pi/2.*sign(relative_pixel(ii,1));
    else
        shift_angle(ii)=atan(relative_pixel(ii,1)./relative_pixel(ii,2));
        if relative_pixel(ii,2)<0
            shift_angle(ii)=shift_angle(ii)-pi;
        end
    end
    
    if r_shift(ii)==0
        f_r(ii)=0;
    else
        f_r(ii)=xsize./r_shift(ii);
    end
    
end

if size_flag==1
    final_xsize=o_xsize;
    final_ysize=o_ysize;
else
    final_xsize=xsize;
    final_ysize=ysize;
end
shift_im=zeros(final_xsize,final_ysize,num);
for ii=1:num
    fr_temp=f_r(ii);
    if fr_temp~=0
        my_angle=shift_angle(ii);
        ft=fftshift(fft2(im(:,:,ii)));
        ft=ft.*exp(-1i*2*pi*(xr.*sin(my_angle)+yr.*cos(my_angle))/fr_temp);
        temp=ifft2(ifftshift(ft));
        shift_im(:,:,ii)=temp(1:final_xsize,1:final_ysize);
    else
        shift_im(:,:,ii)=im(1:final_xsize,1:final_ysize,ii);
    end
end


end

