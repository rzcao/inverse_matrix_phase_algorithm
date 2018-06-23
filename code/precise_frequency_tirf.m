function [ precise_shift,test_ini_phase ] = precise_frequency_tirf(noiseimagef,shiftvalue,search_range)
%estimate the modulation vector k0 down to subpixel accuracy
% return value: precise_shift: the vector k0;
% test_ini_phase: the initial phase.

search_angle_num=10;% search k0 from 10 different angles.
threshold=0.03;% desired accuracy

precise_shift=zeros(size(shiftvalue));
[xsize,ysize,a_num,p_num]=size(noiseimagef);
myangle=0:2*pi/search_angle_num:(2*pi-2*pi/search_angle_num);
sin_angle=sin(myangle);
cos_angle=cos(myangle);
% n_factor=0.0008;% regularization number
test_ini_phase=zeros(1,a_num);
[Y,X]=meshgrid(1:ysize,1:xsize);
xc=floor(xsize/2+1);% center of x axis
yc=floor(ysize/2+1);% center of y axis
yr=Y-yc;
xr=X-xc;
G = 1 - exp(-0.09*sqrt(xr.^2+yr.^2).^1.2);


re_f=zeros(xsize,ysize,a_num,3);
%% 
for ii=1:a_num
    modulation_matrix=[1,1/2*exp(-1i*(pi*0)),1/2*exp(1i*(pi*0));...
                       1,1/2*exp(-1i*(pi*2/3)),1/2*exp(1i*(pi*2/3));...
                       1,1/2*exp(-1i*(pi*4/3)),1/2*exp(1i*(pi*4/3))];
    finv=inv(modulation_matrix);
        re0_temp=zeros(xsize,ysize);
        rep_temp=zeros(xsize,ysize);
%         rem_temp=zeros(xsize,ysize);
        for jj=1:p_num
            re0_temp=finv(1,jj)*noiseimagef(:,:,ii,jj)+re0_temp;
            rep_temp=finv(2,jj)*noiseimagef(:,:,ii,jj)+rep_temp;
%             rem_temp=finv(3,jj)*noiseimagef(:,:,ii,jj)+rem_temp;
        end
        re_f(:,:,ii,1)=re0_temp;
        re_f(:,:,ii,2)=rep_temp;
%         re_f(:,:,ii,3)=rem_temp;
end

norm_ft=zeros(size(re_f));
for ii=1:a_num
    for jj=1:2
        ft_max=max(max(abs(re_f(:,:,ii,jj))));
        re_f(:,:,ii,jj)=re_f(:,:,ii,jj)/ft_max;
%         norm_ft(:,:,ii,jj)=re_f(:,:,ii,jj)./(n_factor+abs(re_f(:,:,ii,jj)));
        norm_ft(:,:,ii,jj)=re_f(:,:,ii,jj).*G;
    end
end

% sub_value=zeros(a_num,3);
% for ii=1:a_num
%     for jj=1:3
%         sub_value(ii,jj)=sum(sum(conj(re_f(:,:,ii,1)).*re_f(:,:,ii,jj)));
%     end
% end

%% local search
for ii=1:a_num
    temp_value=zeros(1,search_angle_num+1);
    temp_phase=zeros(1,search_angle_num+1);
    temp_search=search_range;
    center_x=-shiftvalue(ii,2,1);
    center_y=-shiftvalue(ii,2,2);%
    mytemp=sum(sum(conj(norm_ft(:,:,ii,1)).*circshift(squeeze(norm_ft(:,:,ii,2)),[center_x,center_y])));
    test_ini_phase(1,ii)=angle(mytemp);
    temp_value_max=mytemp;
    temp_phase(1,search_angle_num+1)=angle(mytemp);
    while temp_search>threshold
        for kk=1:search_angle_num
            temp=exact_shift(squeeze(norm_ft(:,:,ii,2)),[center_x+cos_angle(kk)*temp_search,...
                center_y+sin_angle(kk)*temp_search],1);
            mytemp=sum(sum(conj(norm_ft(:,:,ii,1)).*temp));%
            temp_phase(1,kk)=angle(mytemp);
            temp_value(1,kk)=mytemp;
        end
        temp_value(1,kk+1)=temp_value_max;
        temp_max=max(abs(temp_value));
        [xx,yy]=find(abs(temp_value)==temp_max);
        my_num=max(xx,yy);
        temp_value_max=temp_value(my_num);
        if my_num<search_angle_num+0.5
            center_x=center_x+cos_angle(my_num)*temp_search;
            center_y=center_y+sin_angle(my_num)*temp_search;
            test_ini_phase(1,ii)=temp_phase(my_num);
            temp_phase(1,kk+1)=temp_phase(my_num);
        end
        temp_search=temp_search*0.5;
    end
    
    precise_shift(ii,2,1)=-center_x;
    precise_shift(ii,2,2)=-center_y;
    precise_shift(ii,3,1)=center_x;
    precise_shift(ii,3,2)=center_y;
end
        
            


end

