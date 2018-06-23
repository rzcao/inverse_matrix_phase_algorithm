function [solve_ans,solve_sin,solve_cos] = solve_trigonometric_linear_equation_var2(parameter_matrix)
%Equation: a1*cos(x)+b1*sin(x)+c1*cos(y)+d1*sin(y)=0;
%          a2*cos(x)+b2*sin(x)+c2*cos(y)+d2*sin(y)=0;
% Corresponding matrix: [a1,b1,c1,d1;a2,b2,c2,d2];
% prerequisite: a1*a2~=0;||c1*c2~=0;
% M_1sub2=M1-M2;
change_flag=0;
if parameter_matrix(1,1)*parameter_matrix(2,1)==0
    temp(:,1:2)=parameter_matrix(:,3:4);
    temp(:,3:4)=parameter_matrix(:,1:2);
    parameter_matrix=temp;
    change_flag=1; 
end
pm(1,:)=parameter_matrix(1,:)./parameter_matrix(1,1);
pm(2,:)=parameter_matrix(2,:)./parameter_matrix(2,1);

M_1sub2=pm(1,:)-pm(2,:);

coeff_sinx=-1/M_1sub2(2).*[M_1sub2(3),M_1sub2(4)];
coeff_cosx=-pm(1,2).*coeff_sinx-[pm(1,3),pm(1,4)];
% sinx=m1*cosy+m2*siny; coeff_sinx=[m1,m2];
% cosy=l1*cosy+l2*siny; coeff_cosx=[l1,l2];(L)

% eq: 0.5*(m1^2+l1^2-m2^2-l2^2)*cos(2y)+0.5*(m1^2+l1^2+m2^2+l2^2)
%      +(m1*m2+l1*l2)*sin(2y)=1   Acos(2y)+B+Csin(2y)=1;
square_coeff_sin=coeff_sinx.^2;
square_coeff_cos=coeff_cosx.^2;
A=0.5*(square_coeff_sin(1)-square_coeff_sin(2)+square_coeff_cos(1)-square_coeff_cos(2));
B=0.5*(square_coeff_sin(1)+square_coeff_sin(2)+square_coeff_cos(1)+square_coeff_cos(2));
C=coeff_sinx(1)*coeff_sinx(2)+coeff_cosx(1)*coeff_cosx(2);

% eq: (A^2+C^2)cos^2(2y)-2A(1-B)cos(2y)+(1-B)^2-C^2=0;

a=(A^2+C^2);
b=-2*A*(1-B);
c=(1-B)^2-C^2;
delta=b^2-4*a*c;

cos2y_ans(1)=1/(2*a)*(-b+sqrt(delta));
cos2y_ans(2)=1/(2*a)*(-b-sqrt(delta));

if delta>0
    s_v=0;
    if abs(cos2y_ans(1))<=1
        a_cos2y(1)=acos(cos2y_ans(1));
        a_cos2y(2)=-a_cos2y(1);
        s_v=2;
    end
    
    if abs(cos2y_ans(2))<=1
        a_cos2y(s_v+1)=acos(cos2y_ans(2));
        a_cos2y(s_v+2)=-a_cos2y(s_v+1);
    end
end

if (abs(cos2y_ans(2))<=1||abs(cos2y_ans(1))<=1)&&(delta>0)
    length_ans=length(a_cos2y);
    ans_cosy_siny=zeros(2,2*length_ans);
    ans_y=zeros(1,2*length_ans);
    % ans_cosy_siny=[cosy1,cosy2,... ; siny1,siny2,... ];
    % ans_y=[y1,y2,...];
    for ii=1:length_ans
        ans_y(1,2*ii-1)=a_cos2y(ii)*0.5;
        ans_y(1,2*ii)=a_cos2y(ii)*0.5+pi;
    end
    ans_cosy_siny(1,:)=cos(ans_y);
    ans_cosy_siny(2,:)=sin(ans_y);

    ans_cosx=coeff_cosx*ans_cosy_siny;
    ans_sinx=coeff_sinx*ans_cosy_siny;
    ans_x=acos(ans_cosx).*sign(ans_sinx);
else
    length_ans=0.5;
    ans_x=pi;
    ans_y=0;
    ans_sinx=0;
    ans_cosx=-1;
    ans_cosy_siny(1,1)=1;
    ans_cosy_siny(2,1)=0;
end

solve_ans=zeros(2,2*length_ans);
solve_sin=solve_ans;
solve_cos=solve_ans;

if change_flag==0
    solve_ans(1,:)=ans_x;
    solve_ans(2,:)=ans_y;
    solve_sin(1,:)=ans_sinx;
    solve_sin(2,:)=ans_cosy_siny(2,:);
    solve_cos(1,:)=ans_cosx;
    solve_cos(2,:)=ans_cosy_siny(1,:);
else
    solve_ans(2,:)=ans_x;
    solve_ans(1,:)=ans_y;
    solve_sin(2,:)=ans_sinx;
    solve_sin(1,:)=ans_cosy_siny(2,:);
    solve_cos(2,:)=ans_cosx;
    solve_cos(1,:)=ans_cosy_siny(2,:);
end

end

