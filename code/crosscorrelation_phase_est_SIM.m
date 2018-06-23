function [ans_phase] = crosscorrelation_phase_est_SIM(noiseimagef,precise_shift,sigma,OTF)
% Optimize the phases with cross-corrlation method
% This cross-correlation based algorithm is proposed by Wicker et al.
% Sigma is the standard deviation of the guassian noise
%
% The function is implemented based on the following paper:
% K. Wicker, O. Mandula, G. Best, R. Fiolka, and R. Heintzmann,
% "Phase optimisation for structured illumination microscopy," Opt. Express 21, 2032-2049 (2013).
%
% minFunc is available here:
% M. Schmidt. minFunc: unconstrained differentiable multivariate optimization in Matlab.
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.
[~,~,a_num,p_num]=size(noiseimagef);
tensor_d=zeros(p_num,p_num);
ans_phase=zeros(a_num,p_num);
tic;
epsilon=0.01;
for ii=1:a_num
    temp=exact_shift(OTF,precise_shift(ii,2,:),1);
    weigh_mask=temp.*conj(OTF)./(sigma^2*(temp.*conj(temp)+OTF.*conj(OTF))+epsilon);
    
    shift_temp=repmat(squeeze(precise_shift(ii,2,:))',p_num,1);
    temp_ft=exact_shift(squeeze(noiseimagef(:,:,ii,:)),shift_temp,1);
    
    for jj=1:p_num
        for kk=jj:p_num
            tensor_d(jj,kk)=...
                sum(sum(weigh_mask.*noiseimagef(:,:,ii,jj).*conj(temp_ft(:,:,kk))))/...
                sum(weigh_mask(:));
        end
        
        if jj-1>0
            for kk=1:jj-1
                tensor_d(jj,kk)=tensor_d(kk,jj);
            end
        end
    end
    
    ini_phase=[0,pi*2/3,pi*4/3]';
    cost_func=@(est)cost_func_CC(tensor_d,est);
    
    options.display = 'none';
    options.maxFunEvals = 5000;
    options.numDiff=1;
    options.Method = 'lbfgs';
    options.progTol=1e-150;
    [temp,~,~,output_value]=minFunc(cost_func,ini_phase,options);
    
    
    modulation_matrix=[1,1/2*exp(-1i*temp(1,1)),1/2*exp(1i*temp(1,1));...
                       1,1/2*exp(-1i*temp(2,1)),1/2*exp(1i*temp(2,1));...
                       1,1/2*exp(-1i*temp(3,1)),1/2*exp(1i*temp(3,1))];
    matrix_inv=inv(modulation_matrix);
    transpose_inv=matrix_inv';
    result_mat=(modulation_matrix\tensor_d)*transpose_inv;
    global_phase=angle(result_mat(1,3));
    
    ans_phase(ii,:)=temp'-global_phase;
end
runtime_CC=toc;
end