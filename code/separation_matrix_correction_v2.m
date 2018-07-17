function [my_phase] = separation_matrix_correction_v2(noiseimagef,precise_shift,OTF)
% solve the function f(phi1,phi2,phi3)=phase_peak
% version 2 is typically slower than version 3. But it is easier for others
% to understand the algorithm using this version.

tic;

xsize,ysize,a_num,p_num]=size(noiseimagef);
[Y,X]=meshgrid(1:ysize,1:xsize);
xc=floor(xsize/2+1);
yc=floor(ysize/2+1);
yr=Y-yc;% center of x axis
xr=X-xc;% center of y axis
n_filt = 1 - exp(-0.01*sqrt(xr.^2+yr.^2).^1.2);
mi=0.5;

sep_set=[0,2/3*pi,4/3*pi;
    0,1/3*pi,2/3*pi;
    0,4/7*pi,10/7*pi;
    0,4/5*pi,8/5*pi];
sep_num=4;

finv=zeros(3,3,sep_num);
phase_result=zeros(a_num,sep_num);
re_f=zeros(xsize,ysize,a_num,2);
f_true=re_f;

mask=abs(OTF)/(abs(OTF)+0.1);
mask=mask.*n_filt;
for ii=1:a_num
    for jj=1:p_num
        noiseimagef(:,:,ii,jj)=noiseimagef(:,:,ii,jj).*mask;
    end
end

% analytical solution
wide_temp=squeeze(sum(noiseimagef,4)/p_num);
wide_temp=squeeze(sum(wide_temp,3)/a_num);
a=[0.1,0.2,0.15,-0.2];
b=[0.1,-0.1,-0.15,0.2];
c=[-0.1,-0.2,0.3,0.3];
% analytical solution

for my_num=1:sep_num
    %% separation
    for ii=1:a_num
        t_phase(1)=sep_set(my_num,1);
        t_phase(2)=sep_set(my_num,2);
        t_phase(3)=sep_set(my_num,3);
        f_d=[1+a(my_num),mi*exp(-1i*(t_phase(1))),mi*exp(1i*(t_phase(1)));...
            1+b(my_num),mi*exp(-1i*(t_phase(2))),mi*exp(1i*(t_phase(2)));...
            1+c(my_num),mi*exp(-1i*(t_phase(3))),mi*exp(1i*(t_phase(3)))];
        finv(:,:,my_num)=inv(f_d);
        re0_temp=zeros(xsize,ysize);
        rep_temp=zeros(xsize,ysize);
%         rem_temp=zeros(xsize,ysize);
        for jj=1:p_num
            re0_temp=finv(1,jj,my_num)*noiseimagef(:,:,ii,jj)+re0_temp;
            rep_temp=finv(2,jj,my_num)*noiseimagef(:,:,ii,jj)+rep_temp;
%             rem_temp=finv(3,jj,my_num)*noiseimagef(:,:,ii,jj)+rem_temp;
        end
        
        sum_2=sum(squeeze(finv(2,:,my_num)));
        rep_temp=rep_temp-sum_2*wide_temp;
        
        re_f(:,:,ii,1)=re0_temp/3;
        re_f(:,:,ii,2)=rep_temp/3;
        
%         re_f(:,:,ii,1)=re0_temp/3.*n_filt;
%         re_f(:,:,ii,2)=rep_temp/3.*n_filt;
%         re_f(:,:,ii,3)=rem_temp/3.*n_filt;
    end
    
    %% reconstruction
    f_true(:,:,:,1)=re_f(:,:,:,1);
    f_true(:,:,:,2)=exact_shift(squeeze(re_f(:,:,:,2)),-squeeze(precise_shift(:,2,:)),1);
    
    f_reference=sum(f_true(:,:,:,1),3);
    
    for ii=1:a_num 
        phase_result(ii,my_num)=sum(sum(conj(f_reference).*f_true(:,:,ii,2)));
    end
    
end

%% inverse matrix based phase estimation algorithm

inv_matrix=zeros(4,3);
analytical_phase_ans=zeros(a_num,3);
for ii=1:4
    inv_matrix(ii,:)=finv(2,:,ii);
end

for ii=1:a_num
    [~,final_ans]=solve_trigonometric_linear_equation_var3(phase_result(ii,:),inv_matrix);
    analytical_phase_ans(ii,:)=real(final_ans');
end
run_time_var3=toc;
%% inverse matrix based algorithm end
my_phase=mod(analytical_phase_ans,2*pi);

end
