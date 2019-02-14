function [ CR ] = confidence_region( p_fit )

exp_data = xlsread('Dataset1');
time = exp_data(:,1);
Vc0 = exp_data(1,2:end);

[t,Vc_fit] = ode15s( @(t,Vc)sensitivity_analysis_ODEs( t,Vc,p_fit ), time, [Vc0 zeros(1,27*8)]);

%%

%Calculate Parameter Covariance


S_all = [];
V_all = [];
for k = 1:length(time)
    
    S = reshape(Vc_fit(k,28:end),27,8);
    S_all = [S_all; S];
    
    V = (exp_data(k,2:end).*0.05).^2;
    if k == 1
        V(:) = min(V(V>0));  
    end
    V = diag(V);
    V_all = mdiag(V_all,V);
    
end

% [U,E,V] = svd(S_all)

Cov_p = inv(S_all'*inv(V_all)*S_all);

%Parameter Confidence Intervals 97.5%

ndata = length(exp_data(2:end,2:end))*length(exp_data(2:end,1));
np = length(p_fit);
sigma = sqrt(diag(Cov_p));

CR_lb = p_fit' - sigma*tinv(0.975,ndata-np);
CR_ub = p_fit' + sigma*tinv(0.975,ndata-np);
CR = [CR_lb CR_ub]

