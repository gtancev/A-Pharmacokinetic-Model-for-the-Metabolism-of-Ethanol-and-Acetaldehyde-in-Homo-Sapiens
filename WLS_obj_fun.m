function PHI = WLS_obj_fun(p,time,exp_data,W)
    
    Vc0 = exp_data(1,2:end);                 % intial condition of all metabolites
    
    [t,Vc_fit] = ode15s( @(t,Vc)model_odes(t,Vc,p), time, Vc0 ); % solution of the ODE system
    
    wls_obj_fun_value = 0;
    
    for metabolite_number = 1:1:27
        
        Vc_exp = [];            
        Vc_exp = exp_data(2:end,(metabolite_number+1)); % first column of Expt Data is the time
        
        wls_obj_fun_value = wls_obj_fun_value + ( W(metabolite_number).*( Vc_exp - Vc_fit(2:end,metabolite_number) )' * (W(metabolite_number).*( Vc_exp - Vc_fit(2:end,metabolite_number) ) ));        
        
    end
      
    PHI = wls_obj_fun_value;

end