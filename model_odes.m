function dVc = model_odes(t,Vc,p)

k_s_max = p(1);             % maximum stomach emtying rate constant
a = p(2);                   % empiric parameter (dependend on the dose)
V_maxAl = p(3);             % 
V_revAl = p(4);
K_Al = p(5);
K_revAl = p(6);
V_maxAc = p(7);
K_Ac = p(8);

%--------------------------------------------------------------------------
% Tissue Water Volumes


V_C = 11.56;                        % [L] of central compartment
V_M = 25.8;                        % [L] of muscle compartment
V_G = 2.4;                          % [L] of gastrointestinal compartment
V_L = 1.1;                         % [L] of liver compartment

%--------------------------------------------------------------------------
% other parameters

N = 10; 
V_s0 = 0.15;                                % [mL]
C_SAl = 1623;%V_proc*rho_EtOH/46*1000;      % [mM]
D = V_s0*C_SAl;                             % [mmol]

%--------------------------------------------------------------------------
% Important Blood Flowrates

v_L = 1.35;                         % [L/min] in liver (total coming from portal vein and artery)
v_M = 1.5;                          % [L/min]

deltaV = V_L/N;

%--------------------------------------------------------------------------
V_s = Vc(1);
C_GAl = Vc(2);
C_GAc = Vc(3);
C_LAl(1:N) = Vc(4:(N+3));
C_LAc(1:N) = Vc((N+4):(2*N+3));
C_CAl = Vc(2*N+4);
C_CAc = Vc(2*N+5);
C_MAl = Vc(2*N+6);
C_MAc = Vc(2*N+7);

%--------------------------------------------------------------------------
% Calculate Current Stomach Emmptying Rate Constant
k_s = stomach_emptying_const(k_s_max,D,a);

%--------------------------------------------------------------------------
% Derivatives

%-------------------------------
% Stomach Compartement
dV_s = -k_s*V_s;                                            % V_s is volume of fluid in the stomach [L]

%-------------------------------
% Gastrointestinal Compartment
V_GdC_GAl = ((2/3*v_L)*(C_CAl - C_GAl) +  k_s*(V_s)*(C_SAl))/V_G; % mass balance for ethanol
V_GdC_GAc = ((2/3*v_L)*(C_CAc - C_GAc))/V_G;                      % mass balance for acetaldehyde

%-------------------------------
% Liver Compartment 1 (L1)
C_LAl0 = 1/3*C_CAl + 2/3*C_GAl;                              % boundary condition ethanol
C_LAc0 = 1/3*C_CAc + 2/3*C_GAc;                              % boundary condition acetaldehyde

deltaVdC_LAl(1) = (v_L*(C_LAl0 - C_LAl(1)) + (-V_maxAl*C_LAl(1)+V_revAl*C_LAc(1))/(K_Al + C_LAl(1) + K_revAl * C_LAc(1))*deltaV)/deltaV;
deltaVdC_LAc(1) = (v_L*(C_LAc0 - C_LAc(1)) - (-V_maxAl*C_LAl(1)+V_revAl*C_LAc(1))/(K_Al + C_LAl(1) + K_revAl * C_LAc(1))*deltaV + (-V_maxAc*C_LAc(1))/(K_Ac + C_LAc(1))*deltaV)/deltaV;

for i=2:N
deltaVdC_LAl(i) = (v_L*(C_LAl(i-1) - C_LAl(i)) + (-V_maxAl*C_LAl(i)+V_revAl*C_LAc(i))/(K_Al + C_LAl(i) + K_revAl * C_LAc(i))*deltaV)/deltaV;
deltaVdC_LAc(i) = (v_L*(C_LAc(i-1) - C_LAc(i)) - (-V_maxAl*C_LAl(i)+V_revAl*C_LAc(i))/(K_Al + C_LAl(i) + K_revAl * C_LAc(i))*deltaV + (-V_maxAc*C_LAc(i))/(K_Ac + C_LAc(i))*deltaV)/deltaV;
end

%-------------------------------
% Central Compartment
V_CdC_CAl = (-v_L*(C_CAl-C_LAl(N)) - v_M*(C_CAl - C_MAl))/V_C;       % mass balance for ethanol
V_CdC_CAc = (-v_L*(C_CAc-C_LAc(N)) - v_M*(C_CAc - C_MAc))/V_C;       % mass balance for wcetaldehyde

%-------------------------------
% Muscle Compartment
V_MdC_MAl = (v_M*(C_CAl - C_MAl))/V_M;                            % mass balance for ethanol
V_MdC_MAc = (v_M*(C_CAc - C_MAc))/V_M;                            % mass balance for wcetaldehyde

dVc = [dV_s; V_GdC_GAl; V_GdC_GAc; deltaVdC_LAl(:); deltaVdC_LAc(:); V_CdC_CAl; V_CdC_CAc; V_MdC_MAl; V_MdC_MAc];
