%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of the Model

%close all; clc; 
clear all;
p(1) = 0.060; % k_s_max
p(2) = 1.5;   % a
p(3) = 2.2;   % v_MAXaL
p(4) = 29.1;    % V_revAl
p(5) = 0.3898;  % K_Al
p(6) = 1;       % K_revAl
p(7) = 2.74;    % V_maxAc
p(8) = 0.0015;  % K_Ac

t_span = [0 180];

N = 10;
Vc0 = [0.15 zeros(1,2*N+6)];

sol = ode15s(@(t,c)model_odes(t,c,p),t_span,Vc0);%,options);

figure
plot(sol.x,sol.y(1,:))
title('Liquid Volume in Stomach')

figure
plot(sol.x,sol.y([2 3],:))
title('Gastricular Tract')
legend('C_{Al}','C_{Ac}')

figure
plot(sol.x,sol.y([N+3 2*N+3],:))
title('Liver')
legend('C_{Al}','C_{Ac}')

figure
plot(sol.x,sol.y([2*N+4],:))
title('Ethanol Conc. in Central Fluid')
xlabel('Time [min]')
ylabel('Conc. [mM]')

figure
plot(sol.x,sol.y([2*N+5],:).*1000)
title('Acetaldehyde Conc. in Central Fluid')
xlabel('Time [min]')
ylabel('Conc. [\muM]')

figure
plot(sol.x,sol.y([2*N+6 2*N+7],:))
title('Muscle')
legend('C_{Al}','C_{Ac}')