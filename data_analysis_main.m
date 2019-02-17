%% code for data analysis (ML)

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

p_initial = [0.5257 71.7950 1.8979 62.7180 0.3850 80.2816 2.7987 0.0016];                                  % random values between 0 and 40 for each parameter
A_ineq = []; b_ineq = []; A_eq = []; b_eq = [];         % no inequality and equality constraints 
lb = 0*ones(total_parameters,1);                        % lower bound
ub = 100*ones(total_parameters,1);                       % upper bound
total_starting_points = 10;

%% optimisation (maximum likelihood)

problem = createOptimProblem('fmincon', 'objective',@(p)ML_obj_fun(p,time,exp_data), 'lb',lb, 'ub',ub, 'x0',p_initial );
ms = MultiStart;

[p3,fval3] = run(ms,problem,total_starting_points);

%% evaluation of the solution

Vc0 = exp_data(1,2:end);
sol1 = ode15s(@(t,c)model_odes(t,c,p3),time,Vc0);
fit_at_time = deval(sol1,exp_data(:,1));
%% figures

% stomach
figure1 = figure('PaperOrientation','landscape');
axes1 = axes('Parent',figure1,...
    'Position',[0.060313630880579 0.137946127946128 0.869722557297949 0.787053872053872]);
plot(sol1.x,sol1.y(1,:),'-',exp_data(:,1),exp_data(:,2),'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
xlabel('time / min');
ylabel('volume / L');
legend({'V_{Al,fit}','V_{Al,exp}'},'Box','off','FontSize',14);

% stomach residuals

figure
stomach_resid = exp_data(:,2)-fit_at_time(1,:)';
plot(fit_at_time(1,:),stomach_resid,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
%legend({'V_{Al,fit}','V_{Al,exp}'},'Box','off','FontSize',14);
title('Residuals Stomach Compartment (Beverage Volume)')


% gastricular tract
figure2 = figure('PaperOrientation','landscape');
axes2 = axes('Parent',figure2,...
    'Position',[0.060313630880579 0.137946127946128 0.869722557297949 0.787053872053872]);
yyaxis left;
plot(sol1.x,sol1.y([2],:),exp_data(:,1),exp_data(:,3),'o','Linewidth',1.5);
ylabel('concentration / mM');
yyaxis right;
plot(sol1.x,sol1.y([3],:).*10^3,exp_data(:,1),exp_data(:,4).*10^3,'o','Linewidth',1.5);
ylabel('concentration / \muM');
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
xlabel('time / min');
legend({'C_{Al,fit}','C_{Al,exp}','C_{Ac,fit}','C_{Ac,exp}'},'Box','off','FontSize',14);

% Residuals calc and plot: Gastrointestinal
gast_resid_Et = exp_data(:,3)-fit_at_time(2,:)';
gast_resid_Ac = exp_data(:,4)-fit_at_time(3,:)';

figure
plot(fit_at_time(2,:),gast_resid_Et,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Gastointestinal Tract Residuals (Ethanol)')

figure
plot(fit_at_time(3,:)',gast_resid_Ac.*10^3,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Gastointestinal Tract Residuals (Acetaldehyde)')



% liver
figure3 = figure('PaperOrientation','landscape');
axes3 = axes('Parent',figure3,...
    'Position',[0.060313630880579 0.137946127946128 0.869722557297949 0.787053872053872]);
yyaxis left;
plot(sol1.x,sol1.y([N+3],:),exp_data(:,1),exp_data(:,14),'o','Linewidth',1.5);
ylabel('concentration / mM');
yyaxis right;
plot(sol1.x,sol1.y([2*N+3],:).*10^3,exp_data(:,1),exp_data(:,24).*10^3,'o','Linewidth',1.5);
ylabel('concentration / \muM');
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
xlabel('time / min');
legend({'C_{Al,fit}','C_{Al,exp}','C_{Ac,fit}','C_{Ac,exp}'},'Box','off','FontSize',14);

% Residuals calc and plot: Liver
liver_resid_Et = exp_data(:,14)-fit_at_time(13,:)';
liver_resid_Ac = exp_data(:,24)-fit_at_time(23,:)';

figure
plot(fit_at_time(13,:),liver_resid_Et,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Liver Tract Residuals (Ethanol)')

figure
plot(fit_at_time(23,:),liver_resid_Ac.*10^3,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Liver Tract Residuals (Acetaldehyde)')



% central fluid 
figure4 = figure('PaperOrientation','landscape');
axes4 = axes('Parent',figure4,...
    'Position',[0.060313630880579 0.137946127946128 0.869722557297949 0.787053872053872]);
yyaxis left;
plot(sol1.x,sol1.y([2*N+4],:),exp_data(:,1),exp_data(:,25),'o','Linewidth',1.5);
ylabel('concentration / mM');
yyaxis right;
plot(sol1.x,sol1.y([2*N+5],:).*10^3,exp_data(:,1),exp_data(:,26).*10^3,'o','Linewidth',1.5);
ylabel('concentration / \muM');
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
xlabel('time / min');
title('concentrations of ethanol & acetaldehyde in the central fluid')
legend({'C_{Al,fit}','C_{Al,exp}','C_{Ac,fit}','C_{Ac,exp}'},'Box','off','FontSize',14);

% Residuals calc and plot: Central Fluid
central_resid_Et = exp_data(:,25)-fit_at_time(24,:)';
central_resid_Ac = exp_data(:,26)-fit_at_time(25,:)';

figure
plot(fit_at_time(24,:),central_resid_Et,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Central Compartment Residuals (Ethanol)')

figure
plot(fit_at_time(25,:),central_resid_Ac.*10^3,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Central Compartment Residuals (Acetaldehyde)')



% muscle
figure5 = figure('PaperOrientation','landscape');
axes5 = axes('Parent',figure5,...
    'Position',[0.060313630880579 0.137946127946128 0.869722557297949 0.787053872053872]);
yyaxis left;
plot(sol1.x,sol1.y([2*N+6],:),exp_data(:,1),exp_data(:,27),'o','Linewidth',1.5);
ylabel('concentration / mM');
yyaxis right;
plot(sol1.x,sol1.y([2*N+7],:).*10^3,exp_data(:,1),exp_data(:,28).*10^3,'o','Linewidth',1.5);
ylabel('concentration / \muM');
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
xlabel('time / min');
legend({'C_{Al,fit}','C_{Al,exp}','C_{Ac,fit}','C_{Ac,exp}'},'Box','off','FontSize',14);

% Residuals calc and plot; Muscle
muscle_resid_Et = exp_data(:,27)-fit_at_time(26,:)';
muscle_resid_Ac = exp_data(:,28)-fit_at_time(27,:)';

figure
plot(fit_at_time(26,:),muscle_resid_Et,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Muscle Compartment Residuals (Ethanol)')

figure
plot(fit_at_time(27,:),muscle_resid_Ac.*10^3,'o','Linewidth',1.5);
set(gca,'FontSize',16,'Linewidth',1.5);
set(gcf,'Position',[10,10,1800,600]);
ylabel('y_m - y');
xlabel('y');
title('Muscle Compartment Residuals (Acetaldehyde)')



%% Identifiability Analysis

%%

S_initial = [Vc0 zeros(1,27*8)];
[t,Vc_fit] = ode15s( @(t,Vc)sensitivity_analysis_ODEs( t,Vc,p3 ), time, S_initial);

%%

% Calculate Parameter Covariance


S_all = [];
V_all = [];
for k = 2:length(time)
    
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
% rank(S);

Cov_p = inv(S_all'*inv(V_all)*S_all);

% Parameter Confidence Intervals 97.5%

ndata = length(exp_data(2:end,2:end))*length(exp_data(2:end,1));
np = length(p3);
sigma = sqrt(abs(diag(Cov_p)));

CR = [p3' - sigma*tinv(0.975,ndata-np), p3' + sigma*tinv(0.975,ndata-np)];