function PHI = OLS_obj_fun(p,time,exp_data)
    
    Vc0 = exp_data(1,2:end);                 % intial condition of all metabolites
    
    [t,Vc] = ode15s( @(t,Vc)model_odes(t,Vc,p), time, Vc0 ); % solution of the ODE system
    
    ols_obj_fun_value = 0;
    
    for metabolite_number = 1:1:27
        
        Vc_exp = [];            
        Vc_exp = exp_data(2:end,(metabolite_number+1)); % first column of Expt Data is the time
        Vc_fit = Vc(2:end,metabolite_number);
        
        ols_obj_fun_value = ols_obj_fun_value + ( ( Vc_exp - Vc_fit )' * ( Vc_exp - Vc_fit ) );        
        
    end
      
    PHI = ols_obj_fun_value;

end