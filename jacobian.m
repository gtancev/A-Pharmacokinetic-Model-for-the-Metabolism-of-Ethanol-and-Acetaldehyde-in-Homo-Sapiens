%% Sensitivity - Jacobian

%% --------------------------------------------------------------------------
% Tissue Water Volumes


V_C = 11.56;                        % [L] of central compartment
V_M = 25.8;                        % [L] of muscle compartment
V_G = 2.4;                          % [L] of gastrointestinal compartment
V_L = 1.1;                         % [L] of liver compartment

%--------------------------------------------------------------------------
% other parameters

N = 10; 
V_s0 = 0.15;                                % [L]
C_SAl = 1623;%V_proc*rho_EtOH/46*1000;      % [mM]
D = V_s0*C_SAl;                             % [mmol]

%--------------------------------------------------------------------------
% Important Blood Flowrates

v_L = 1.35;                         % [L/min] in liver (total coming from portal vein and artery)
v_M = 1.5;                          % [L/min]

deltaV = V_L/N;

%%

syms k_s_max a V_maxAl V_revAl K_Al K_revAl V_maxAc K_Ac
syms V_s C_GAl C_GAc C_LAl1 C_LAl2 C_LAl3 C_LAl4 C_LAl5 C_LAl6 C_LAl7 C_LAl8 C_LAl9 C_LAl10 C_LAc1 C_LAc2 C_LAc3 C_LAc4 C_LAc5 C_LAc6 C_LAc7 C_LAc8 C_LAc9 C_LAc10 C_CAl C_CAc C_MAl C_MAc

%% 
% Derivatives
%-------------------------------
% Stomach Compartement

k_s = k_s_max/(1 + a*(D/1000)^2);
dV_s = -k_s*V_s;                                            % V_s is volume of fluid in the stomach [L]

%-------------------------------
% Gastrointestinal Compartment

V_GdC_GAl = ((2/3*v_L)*(C_CAl - C_GAl) +  k_s*(V_s)*(C_SAl))/V_G; % mass balance for ethanol
V_GdC_GAc = ((2/3*v_L)*(C_CAc - C_GAc))/V_G;                      % mass balance for acetaldehyde

%-------------------------------
% Liver Compartments 1-10

C_LAl0 = 1/3*C_CAl + 2/3*C_GAl;                              % boundary condition ethanol
C_LAc0 = 1/3*C_CAc + 2/3*C_GAc;                              % boundary condition acetaldehyde

deltaVdC_LAl(1) = (v_L*(C_LAl0 - C_LAl1) + (-V_maxAl*C_LAl1+V_revAl*C_LAc1)/(K_Al + C_LAl1 + K_revAl * C_LAc1)*deltaV)/deltaV;
deltaVdC_LAc(1) = (v_L*(C_LAc0 - C_LAc1) - (-V_maxAl*C_LAl1+V_revAl*C_LAc1)/(K_Al + C_LAl1 + K_revAl * C_LAc1)*deltaV + (-V_maxAc*C_LAc1)/(K_Ac + C_LAc1)*deltaV)/deltaV;

deltaVdC_LAl(2) = (v_L*(C_LAl1 - C_LAl2) + (-V_maxAl*C_LAl2+V_revAl*C_LAc2)/(K_Al + C_LAl2 + K_revAl * C_LAc2)*deltaV)/deltaV;
deltaVdC_LAc(2) = (v_L*(C_LAc1 - C_LAc2) - (-V_maxAl*C_LAl2+V_revAl*C_LAc2)/(K_Al + C_LAl2 + K_revAl * C_LAc2)*deltaV + (-V_maxAc*C_LAc2)/(K_Ac + C_LAc2)*deltaV)/deltaV;

deltaVdC_LAl(3) = (v_L*(C_LAl2 - C_LAl3) + (-V_maxAl*C_LAl3+V_revAl*C_LAc3)/(K_Al + C_LAl3 + K_revAl * C_LAc3)*deltaV)/deltaV;
deltaVdC_LAc(3) = (v_L*(C_LAc2 - C_LAc3) - (-V_maxAl*C_LAl3+V_revAl*C_LAc3)/(K_Al + C_LAl3 + K_revAl * C_LAc3)*deltaV + (-V_maxAc*C_LAc3)/(K_Ac + C_LAc3)*deltaV)/deltaV;

deltaVdC_LAl(4) = (v_L*(C_LAl3 - C_LAl4) + (-V_maxAl*C_LAl4+V_revAl*C_LAc4)/(K_Al + C_LAl4 + K_revAl * C_LAc4)*deltaV)/deltaV;
deltaVdC_LAc(4) = (v_L*(C_LAc3 - C_LAc4) - (-V_maxAl*C_LAl4+V_revAl*C_LAc4)/(K_Al + C_LAl4 + K_revAl * C_LAc4)*deltaV + (-V_maxAc*C_LAc4)/(K_Ac + C_LAc4)*deltaV)/deltaV;

deltaVdC_LAl(5) = (v_L*(C_LAl4 - C_LAl5) + (-V_maxAl*C_LAl5+V_revAl*C_LAc5)/(K_Al + C_LAl5 + K_revAl * C_LAc5)*deltaV)/deltaV;
deltaVdC_LAc(5) = (v_L*(C_LAc4 - C_LAc5) - (-V_maxAl*C_LAl5+V_revAl*C_LAc5)/(K_Al + C_LAl5 + K_revAl * C_LAc5)*deltaV + (-V_maxAc*C_LAc5)/(K_Ac + C_LAc5)*deltaV)/deltaV;

deltaVdC_LAl(6) = (v_L*(C_LAl5 - C_LAl6) + (-V_maxAl*C_LAl6+V_revAl*C_LAc6)/(K_Al + C_LAl6 + K_revAl * C_LAc6)*deltaV)/deltaV;
deltaVdC_LAc(6) = (v_L*(C_LAc5 - C_LAc6) - (-V_maxAl*C_LAl6+V_revAl*C_LAc6)/(K_Al + C_LAl6 + K_revAl * C_LAc6)*deltaV + (-V_maxAc*C_LAc6)/(K_Ac + C_LAc6)*deltaV)/deltaV;

deltaVdC_LAl(7) = (v_L*(C_LAl6 - C_LAl7) + (-V_maxAl*C_LAl7+V_revAl*C_LAc7)/(K_Al + C_LAl7 + K_revAl * C_LAc7)*deltaV)/deltaV;
deltaVdC_LAc(7) = (v_L*(C_LAc6 - C_LAc7) - (-V_maxAl*C_LAl7+V_revAl*C_LAc7)/(K_Al + C_LAl7 + K_revAl * C_LAc7)*deltaV + (-V_maxAc*C_LAc7)/(K_Ac + C_LAc7)*deltaV)/deltaV;

deltaVdC_LAl(8) = (v_L*(C_LAl7 - C_LAl8) + (-V_maxAl*C_LAl8+V_revAl*C_LAc8)/(K_Al + C_LAl8 + K_revAl * C_LAc8)*deltaV)/deltaV;
deltaVdC_LAc(8) = (v_L*(C_LAc7 - C_LAc8) - (-V_maxAl*C_LAl8+V_revAl*C_LAc8)/(K_Al + C_LAl8 + K_revAl * C_LAc8)*deltaV + (-V_maxAc*C_LAc8)/(K_Ac + C_LAc8)*deltaV)/deltaV;

deltaVdC_LAl(9) = (v_L*(C_LAl8 - C_LAl9) + (-V_maxAl*C_LAl9+V_revAl*C_LAc9)/(K_Al + C_LAl9 + K_revAl * C_LAc9)*deltaV)/deltaV;
deltaVdC_LAc(9) = (v_L*(C_LAc8 - C_LAc9) - (-V_maxAl*C_LAl9+V_revAl*C_LAc9)/(K_Al + C_LAl9 + K_revAl * C_LAc9)*deltaV + (-V_maxAc*C_LAc9)/(K_Ac + C_LAc9)*deltaV)/deltaV;

deltaVdC_LAl(10) = (v_L*(C_LAl9 - C_LAl10) + (-V_maxAl*C_LAl10+V_revAl*C_LAc10)/(K_Al + C_LAl10 + K_revAl * C_LAc10)*deltaV)/deltaV;
deltaVdC_LAc(10) = (v_L*(C_LAc9 - C_LAc10) - (-V_maxAl*C_LAl10+V_revAl*C_LAc10)/(K_Al + C_LAl10 + K_revAl * C_LAc10)*deltaV + (-V_maxAc*C_LAc10)/(K_Ac + C_LAc10)*deltaV)/deltaV;

%-------------------------------
% Central Compartment
V_CdC_CAl = (-v_L*(C_CAl-C_LAl10) - v_M*(C_CAl - C_MAl))/V_C;       % mass balance for ethanol
V_CdC_CAc = (-v_L*(C_CAc-C_LAc10) - v_M*(C_CAc - C_MAc))/V_C;       % mass balance for acetaldehyde

%-------------------------------
% Muscle Compartment
V_MdC_MAl = (v_M*(C_CAl - C_MAl))/V_M;                            % mass balance for ethanol
V_MdC_MAc = (v_M*(C_CAc - C_MAc))/V_M;                            % mass balance for wcetaldehyde

%%

dVc = [dV_s; V_GdC_GAl; V_GdC_GAc; deltaVdC_LAl(:); deltaVdC_LAc(:); V_CdC_CAl; V_CdC_CAc; V_MdC_MAl; V_MdC_MAc];

%%

J = jacobian(dVc,[V_s C_GAl C_GAc C_LAl1 C_LAl2 C_LAl3 C_LAl4 C_LAl5 C_LAl6 C_LAl7 C_LAl8 C_LAl9 C_LAl10 C_LAc1 C_LAc2 C_LAc3 C_LAc4 C_LAc5 C_LAc6 C_LAc7 C_LAc8 C_LAc9 C_LAc10 C_CAl C_CAc C_MAl C_MAc])
dfdp = jacobian(dVc,[k_s_max a V_maxAl V_revAl K_Al K_revAl V_maxAc K_Ac])
