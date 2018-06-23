function [my_phase] = separation_matrix_correction_v3(noiseimagef,precise_shift,OTF)
% solve the function f(phi1,phi2,phi3)=phase_peak
% This function calculate the tensor that is similar to the one that proposed
% by Wicker et al [1] to increase the computing speed. See Eq. 11 in Ref [1].
% Here, we use the conventional correlation instead the weighted correaltion.
% And we use the conjugation of the Fourier transform of the estimated
% widefield result in our correlation.
%
% Reference:
% [1] K. Wicker, O. Mandula, G. Best, R. Fiolka, and R. Heintzmann,
% "Phase optimisation for structured illumination microscopy," Opt. Express 21, 2032-2049 (2013).

[~,~,a_num,p_num]=size(noiseimagef);
% [Y,X]=meshgrid(1:ysize,1:xsize);
% xc=floor(xsize/2+1);
% yc=floor(ysize/2+1);
% yr=Y-yc;% center of x axis
% xr=X-xc;% center of y axis
% n_filt = 1 - exp(-0.01*sqrt(xr.^2+yr.^2).^1.2);
mi=0.5;

sep_set=[0,2/3*pi,4/3*pi;
    0,1/3*pi,2/3*pi;
    0,4/7*pi,10/7*pi;
    0,4/5*pi,8/5*pi];
sep_num=4;

redundant=zeros(a_num,1);
ref_wide=zeros(a_num,p_num);
finv=zeros(3,3,sep_num);
phase_result=zeros(a_num,sep_num);

tic;

mask=abs(OTF./(OTF+0.1));
% for ii=1:a_num
%     for jj=1:p_num
%         noiseimagef(:,:,ii,jj)=noiseimagef(:,:,ii,jj).*mask;
%     end
% end

% analytical solution
wide_temp=squeeze(sum(noiseimagef,4)/p_num);
wide_temp=squeeze(sum(wide_temp,3)/a_num);
wide_temp=conj(wide_temp);
a=[0.1,0.2,0.15,-0.2];
b=[0.1,-0.1,-0.15,0.2];
c=[-0.1,-0.2,0.3,0.3];
% analytical solution
for ii=1:a_num
    
    shift_temp=repmat(squeeze(precise_shift(ii,2,:))',p_num,1);
    temp_ft=exact_shift(squeeze(noiseimagef(:,:,ii,:)),-shift_temp,1);
    shift_wide=exact_shift(wide_temp,-shift_temp(1,:),1);
    
    redundant(ii,1)=sum(sum(wide_temp.*conj(shift_wide)));
    for jj=1:p_num
        ref_wide(ii,jj)=sum(sum(wide_temp.*temp_ft(:,:,jj)));
    end
end


for my_num=1:sep_num
    %% separation
    for ii=1:a_num
        t_phase=sep_set(my_num,:);
        
        f_d=[1+a(my_num),mi*exp(-1i*(t_phase(1))),mi*exp(1i*(t_phase(1)));...
            1+b(my_num),mi*exp(-1i*(t_phase(2))),mi*exp(1i*(t_phase(2)));...
            1+c(my_num),mi*exp(-1i*(t_phase(3))),mi*exp(1i*(t_phase(3)))];
        finv(:,:,my_num)=inv(f_d);
        transpose_inv=(finv(:,:,my_num))';
        
        sum_2=sum(squeeze(finv(2,:,my_num)));
        temp_mat=ref_wide(ii,:)*conj(transpose_inv);
        phase_result(ii,my_num)=temp_mat(1,2)-sum_2*redundant(ii,1);
%         temp_mat=(f_d\tensor_d(:,:,ii))*transpose_inv;
%         phase_result(ii,my_num)=temp_mat(1,2)-res_wide(ii,:)*transpose_inv(:,1)*sum_2;
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