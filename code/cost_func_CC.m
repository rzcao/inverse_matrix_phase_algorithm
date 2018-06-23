function [cost] = cost_func_CC(tensor_d,cal_phase)
% The cost function in the cross-correlation based phase estimation method
% 
modulation_matrix=[1,1/2*exp(-1i*cal_phase(1,1)),1/2*exp(1i*cal_phase(1,1));...
                   1,1/2*exp(-1i*cal_phase(2,1)),1/2*exp(1i*cal_phase(2,1));...
                   1,1/2*exp(-1i*cal_phase(3,1)),1/2*exp(1i*cal_phase(3,1))];
matrix_inv=inv(modulation_matrix);
conj_transpose_inv=matrix_inv';

result_mat=(modulation_matrix\tensor_d)*conj_transpose_inv;
abs_mat=abs(result_mat);
% abs_mat=result_mat.*conj(result_mat);
cost=abs_mat(1,1)+abs_mat(2,2)+abs_mat(3,3)...
    +abs_mat(1,2)+abs_mat(3,1);
cost=sqrt(cost);
end

