function [logLik]=Likelihood_SVAR(teta)

global Sigma_u
% global log_lk_reduced
global T
global ParamNumberB

B = zeros(size(Sigma_u,1),size(Sigma_u,1));

    for c_par = 1 : size(ParamNumberB,1)
    B(ParamNumberB(c_par,1))=teta(c_par);   
    end

    K = (B)^(-1);   
    M=size(K,1);

    logLik =-(-0.5*T*M*(log(2*pi))...                               the sign has been changed because the procedure minimizes 
	             + 0.5*T*log((det(K)^2))-0.5*T*trace((K'*K)*Sigma_u));
    

end