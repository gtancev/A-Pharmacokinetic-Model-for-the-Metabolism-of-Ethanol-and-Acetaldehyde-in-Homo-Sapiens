%% code for data analysis (OLS)

clear all; close all; clc;

%% load data

exp_data = xlsread('Dataset1');

%% preparation

time = exp_data(:,1);
total_time_points = length(time);
total_metabolites = length(exp_data(1,:))-1;
total_parameters = 8;
N = 10;

%% configurations for optimisation

p_initial = [0.5257 71.7950 1.8979 62.7180 0.3850 80.2816 2.7987 0.0016];                               % random values between 0 and 40 for each parameter
A_ineq = []; b_ineq = []; A_eq = []; b_eq = [];         % no inequality and equality constraints 
lb = 0*ones(total_parameters,1);                        % lower bound
ub = 100*ones(total_parameters,1);                       % upper bound
total_starting_points = 10;

%% optimisation (weighted least squares, WLS)

W = ones(length(exp_data(1,:))-1,1);
W(3) = sum(exp_data(2:end,3)./exp_data(2:end,4))/length(exp_data(2:end,3));
W(14:23) = sum(exp_data(2:end,13)./exp_data(2:end,23))/length(exp_data(2:end,13));
W(25) = sum(exp_data(2:end,25)./exp_data(2:end,25))/length(exp_data(2:end,25));
W(27) = sum(exp_data(2:end,27)./exp_data(2:end,28))/length(exp_data(2:end,27));
      
problem = createOptimProblem('fmincon', 'objective',@(p)WLS_obj_fun(p,time,exp_data,W), 'lb',lb, 'ub',ub, 'x0',p_initial );
ms = MultiStart;
[p2,fval2] = run(ms,problem,total_starting_points);

%% evaluation of the solution

Vc0 = exp_data(1,2:end);
sol2 = ode15s(@(t,c)model_odes(t,c,p2),time,Vc0);

%% figures

figure(6);
plot(sol2.x,sol2.y(1,:),'-',exp_data(:,1),exp_data(:,2),'o');
set(gca,'FontSize',12);
title('liquid volume in stomach');
xlabel('time / [min]');
ylabel('volume / [L]');

figure(7);
subplot(1,2,1);
plot(sol2.x,sol2.y([2],:),exp_data(:,1),exp_data(:,3),'o');
set(gca,'FontSize',12);
title('gastricular tract');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol2.x,sol2.y([3],:),exp_data(:,1),exp_data(:,4),'o');
set(gca,'FontSize',12);
title('gastricular tract');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

figure(8);
subplot(1,2,1);
plot(sol2.x,sol2.y([N+3],:),exp_data(:,1),exp_data(:,14),'o');
set(gca,'FontSize',12);
title('liver');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol2.x,sol2.y([2*N+3],:),exp_data(:,1),exp_data(:,24),'o');
set(gca,'FontSize',12);
title('liver');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

figure(9)
subplot(1,2,1);
plot(sol2.x,sol2.y([2*N+4],:),exp_data(:,1),exp_data(:,25),'o');
set(gca,'FontSize',12);
title('central fluid')
xlabel('time / [min]')
ylabel('concentration / [mM]')
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol2.x,sol2.y([2*N+5],:),exp_data(:,1),exp_data(:,26),'o');
set(gca,'FontSize',12);
title('central fluid');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

figure(10)
subplot(1,2,1);
plot(sol2.x,sol2.y([2*N+6],:),exp_data(:,1),exp_data(:,27),'o');
set(gca,'FontSize',12);
title('muscle');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol2.x,sol2.y([2*N+7],:),exp_data(:,1),exp_data(:,28),'o');
set(gca,'FontSize',12);
title('muscle');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

%%

% for k=1:total_time_points;
%     W2(:,k) = W;
% end
% 
% W2 = W2';
% 
% %%
% 
% CR = confidence_region( p2 );