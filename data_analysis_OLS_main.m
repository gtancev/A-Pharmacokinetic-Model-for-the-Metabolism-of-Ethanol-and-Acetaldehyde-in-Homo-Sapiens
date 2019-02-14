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

p_initial = 100*rand(8,1);                               % random values between 0 and 40 for each parameter
A_ineq = []; b_ineq = []; A_eq = []; b_eq = [];         % no inequality and equality constraints 
lb = 0*ones(total_parameters,1);                        % lower bound
ub = 100*ones(total_parameters,1);                       % upper bound
total_starting_points = 10;

%% normalisation of data

% for i=1:total_metabolites+1
%    exp_data(:,i) =  exp_data(:,i)/max(exp_data(:,i));
% end

%% optimisation (ordinary least squares, OLS)

problem = createOptimProblem('fmincon', 'objective',@(p)OLS_obj_fun(p,time,exp_data), 'lb',lb, 'ub',ub, 'x0',p_initial );
ms = MultiStart;

[p1,fval1] = run(ms,problem,total_starting_points);

%% evaluation of the solution

Vc0 = exp_data(1,2:end);
sol1 = ode15s(@(t,c)model_odes(t,c,p1),time,Vc0);

%% figures

figure(1);
plot(sol1.x,sol1.y(1,:),'-',exp_data(:,1),exp_data(:,2),'o');
set(gca,'FontSize',12);
title('liquid volume in stomach');
xlabel('time / [min]');
ylabel('volume / [L]');

figure(2);
subplot(1,2,1);
plot(sol1.x,sol1.y([2],:),exp_data(:,1),exp_data(:,3),'o');
set(gca,'FontSize',12);
title('gastricular tract');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol1.x,sol1.y([3],:),exp_data(:,1),exp_data(:,4),'o');
set(gca,'FontSize',12);
title('gastricular tract');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

figure(3);
subplot(1,2,1);
plot(sol1.x,sol1.y([N+3],:),exp_data(:,1),exp_data(:,14),'o');
set(gca,'FontSize',12);
title('liver');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol1.x,sol1.y([2*N+3],:),exp_data(:,1),exp_data(:,24),'o');
set(gca,'FontSize',12);
title('liver');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

figure(4)
subplot(1,2,1);
plot(sol1.x,sol1.y([2*N+4],:),exp_data(:,1),exp_data(:,25),'o');
set(gca,'FontSize',12);
title('central fluid')
xlabel('time / [min]')
ylabel('concentration / [mM]')
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol1.x,sol1.y([2*N+5],:),exp_data(:,1),exp_data(:,26),'o');
set(gca,'FontSize',12);
title('central fluid');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

figure(5)
subplot(1,2,1);
plot(sol1.x,sol1.y([2*N+6],:),exp_data(:,1),exp_data(:,27),'o');
set(gca,'FontSize',12);
title('muscle');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Al,fit}','C_{Al,exp}');

subplot(1,2,2);
plot(sol1.x,sol1.y([2*N+7],:),exp_data(:,1),exp_data(:,28),'o');
set(gca,'FontSize',12);
title('muscle');
xlabel('time / [min]');
ylabel('concentration / [mM]');
legend('C_{Ac,fit}','C_{Ac,exp}');

%%

CR = confidence_region( p1 );