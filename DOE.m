%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of Experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

p = importdata('parameters.mat');
V_all = importdata('V_all.mat');
exp_data = xlsread('Dataset1');

init_cond = zeros(27*9,1);

n_points = 20;
V_min = 0.05;
V_max = 0.75;

%%

vol = linspace(V_min,V_max,n_points);
time = exp_data(:,1);

for i=1:n_points
    
init_cond(1) = vol(i);  % change initial condition of alcoholic beverage in 
                        % stomach to desired to next value
[t_sim,Vc_sim] = ode15s(@(t,Vc)sensitivity_analysis_ODEs(t,Vc,p),time,init_cond);

S_all = [];

for k = 1:length(time)
    
    S = reshape(Vc_sim(k,28:end),27,8);
    S_all = [S_all; S];
    
end

% [U,E,V] = svd(S_all)
% rank(S);

FIM = S_all'*inv(V_all)*S_all; % Fisher Information Matrix

tra_FIM(i) = trace(FIM); % T-optimal
% det_FIM(i) = det(FIM); % D-optimal
end

%% Plot(s)

figure1 = figure('PaperOrientation','landscape');
axes1 = axes('Parent',figure1,...
    'Position',[0.060313630880579 0.137946127946128 0.869722557297949 0.787053872053872]);
plot(vol,tra_FIM,'LineWidth',1.5)
xlabel('volume / L')
ylabel('tr(FIM)')
set(gca,'FontSize',16,'LineWidth',1.5)
set(gcf,'Position',[50 50 1800 600])

% figure(2);
% plot(vol,det_FIM,'LineWidth',1.5)
% xlabel('volume / L')
% ylabel('det(FIM)')
% set(gca,'FontSize',16,'LineWidth',1.5)
% set(gcf,'Position',[50 50 1600 800])
    
