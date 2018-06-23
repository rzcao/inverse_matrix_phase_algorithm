function [ solve_ans,final_ans ] = solve_trigonometric_linear_equation_var3(est_phase,inv_matrix)
% solve trigonometric linear equation with additional inverse matrix
% theta_sub=phase_cross-phase_inv (theta1=phi1-b1);
% theta_sub=[4,:] 4 denote 4 equations
% phase_inv=[b1,b2,b3;d1,d2,d3;...]
% amplitude=[a1,a2,a3;c1,c2,c3;...]
ini_inv_matrix=inv_matrix;
est_phase=est_phase./abs(est_phase);
cal_num_def=2;
solve_ans=zeros(3,8,cal_num_def);
% 2 calculation
for change_num=1:cal_num_def
% 2 calculation

    phase_cross=angle(est_phase);
    amplitude=abs(inv_matrix);

    phase_inv=angle(inv_matrix);
    eq_num=4;

    theta_sub=zeros(eq_num,3);
    for ii=1:eq_num
    theta_sub(ii,:)=phase_cross(ii)-phase_inv(ii,:);
    end
    cos_theta=cos(theta_sub);
    sin_theta=sin(theta_sub);
    % eq:a1*sin(phi1-b1)*cos(x)-a2*cos(phi1-b1)*sin(x)+...=0
    %    c1*sin(phi2-d1)*cos(x)-c2*cos(phi2-d1)*sin(x)+...=0
    % A=[a1*sin(phi1-b1),-a2*cos(phi1-b1),...  ;  c1*sin(phi2-d1),-c2*cos(phi2-d1),...]
    A=zeros(eq_num,3*2);
    for ii=1:eq_num
        for jj=1:3
            A(ii,2*jj-1)=amplitude(ii,jj)*sin_theta(ii,jj);
            A(ii,2*jj)=-amplitude(ii,jj)*cos_theta(ii,jj);
        end
    end

    m1=A(2,1)/A(1,1);
    % eq: sinx1=I32cosx2-I42sinx2+I52cosx3-I62sinx3;
    % coeff_sinx1=[I32,-I42,I52,-I62];
    coeff_sinx1=-(m1*A(1,3:6)-A(2,3:6))./(m1*A(1,2)-A(2,2));
    coeff_cosx1=-(A(1,2)*coeff_sinx1+A(1,3:6))/A(1,1);

    my_matrix=zeros(eq_num-2,2*2);
    my_matrix(1,:)=A(3,1)*coeff_cosx1+A(3,2)*coeff_sinx1+A(3,3:6);
    my_matrix(2,:)=A(4,1)*coeff_cosx1+A(4,2)*coeff_sinx1+A(4,3:6);

    [solve_ans_temp,solve_sin,solve_cos]=solve_trigonometric_linear_equation_var2(my_matrix);
    [~,mylength]=size(solve_sin);
    mysolve=zeros(2*3,mylength);
    % mysclve: odd raw: cos      even raw: sin
    mysolve(3,:)=solve_cos(1,:);
    mysolve(5,:)=solve_cos(2,:);
    mysolve(4,:)=solve_sin(1,:);
    mysolve(6,:)=solve_sin(2,:);

    mysolve(1,:)=coeff_cosx1*mysolve(3:6,:);
    mysolve(2,:)=coeff_sinx1*mysolve(3:6,:);

    abs_solve12=sqrt(mysolve(1,:).^2+mysolve(2,:).^2);
    mysolve(1,:)=mysolve(1,:)./abs_solve12;
    mysolve(2,:)=mysolve(2,:)./abs_solve12;

%     solve_ans(:,:)=zeros(3,mylength);
%     solve_ans(1,:)=acos(mysolve(1,:)).*sign(mysolve(2,:));
%     solve_ans(2:3,:)=solve_ans_temp;

% 2 calculation
    if change_num==1
%         solve_ans(:,1:mylength,change_num)=zeros(3,mylength);
        solve_ans(1,1:mylength,change_num)=acos(mysolve(1,:)).*sign(mysolve(2,:));
        solve_ans(2:3,1:mylength,change_num)=solve_ans_temp;
    else
%         solve_ans(:,1:mylength,change_num)=zeros(3,mylength);
        solve_ans(2,1:mylength,change_num)=acos(mysolve(1,:)).*sign(mysolve(2,:));
        solve_ans(1,1:mylength,change_num)=solve_ans_temp(1,:);
        solve_ans(3,1:mylength,change_num)=solve_ans_temp(2,:);
    end
    inv_matrix_temp=inv_matrix;
    inv_matrix(:,2)=inv_matrix_temp(:,1);
    inv_matrix(:,1)=inv_matrix_temp(:,2);
    
end
% 2 calculation

%% find proper answer
[~,ans_num,cal_num]=size(solve_ans);

%% /1.calculation and average, 2.find optimized answer.

% version1: average
% if cal_num>1
%     cal_num=cal_num+1;
%     solve_ans(:,:,cal_num)=sum(solve_ans,3)/(cal_num-1);
% end

%version2: average with propose
% if cal_num>1
%     cal_num=cal_num+1;
%     solve_ans(1,:,3)=solve_ans(1,:,2);
%     solve_ans(2,:,3)=0.5*(solve_ans(2,:,1)+solve_ans(2,:,2));
%     solve_ans(3,:,3)=solve_ans(3,:,1);
% end
% 
% 
% result_ans=zeros(eq_num,ans_num,cal_num);
% error_temp=result_ans;
% for ii=1:cal_num
% result_ans(:,:,ii)=ini_inv_matrix*exp(1i*solve_ans(:,:,ii));
% end
% result_ans=result_ans./abs(result_ans);
% 
% for ii=1:eq_num
%     error_temp(ii,:,:)=result_ans(ii,:,:)-est_phase(ii);
% end
% 
% error_temp=error_temp.*conj(error_temp);
% 
% error=squeeze(sum(error_temp,1));
% [x_num,y_num]=find(error==min(error(:)));
% 
% final_ans=solve_ans(:,x_num,y_num);

%% 1.calculation and average, 2.find optimized answer.  end/

%% /1.calculation and find optimized answer. 2.average

result_ans=zeros(eq_num,ans_num,cal_num);
error_temp=result_ans;
for ii=1:cal_num
result_ans(:,:,ii)=ini_inv_matrix*exp(1i*solve_ans(:,:,ii));
end
result_ans=result_ans./abs(result_ans);
for ii=1:eq_num
    error_temp(ii,:,:)=result_ans(ii,:,:)-est_phase(ii);
end
error_temp=error_temp.*conj(error_temp);
error=squeeze(sum(error_temp,1));

phase_ans=zeros(3,cal_num+1);
temp_error=zeros(1,cal_num+1);
for ii=1:cal_num
    temp=error(:,ii);
    temp_error(1,ii)=min(error(:,ii));
    [x_num,~]=find(temp==temp_error(1,ii));
    phase_ans(:,ii)=solve_ans(:,x_num(1),ii);
end

for ii=1:3
    if abs(phase_ans(ii,1)-phase_ans(ii,2))<pi
        phase_ans(ii,cal_num+1)=(phase_ans(ii,1)+phase_ans(ii,2))/cal_num;
    else
        phase_ans(ii,cal_num+1)=(phase_ans(ii,1)+phase_ans(ii,2)+2*pi)/cal_num;
    end
end

temp=ini_inv_matrix*exp(1i*phase_ans(:,cal_num+1));
temp=temp./abs(temp)-conj(est_phase');
temp=temp.*conj(temp);
temp_error(1,cal_num+1)=sum(temp);

[~,y_num]=find(temp_error==min(temp_error(:)));
final_ans=phase_ans(:,y_num);
%% 1.calculation and find optimized answer. 2.average  end/

end

