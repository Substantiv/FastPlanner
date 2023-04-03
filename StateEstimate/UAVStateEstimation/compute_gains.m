%% Roll Attitude Hold
% Compute roll attitude hold loop gains: kp_phi, kd_phi, ki_phi
delta_a_max = pi/4;
Va = sqrt(x_trim(4)^2+x_trim(5)^2+x_trim(6)); % m/s

%%%DESIGN PARAMETERS%%%
zeta_phi = 1.5;
e_phi_max = 120*pi/180;
%%%END%%%%%%%%%%%%%%%%%

a_phi1 = -P.rho*Va*P.S_wing*P.b^2*P.C_p_p/4;
a_phi2 = P.rho*Va^2*P.S_wing*P.b*P.C_p_delta_a/2;

wn_phi = sqrt(abs(a_phi2)*delta_a_max/e_phi_max);

P.kp_phi = 2.87; %delta_a_max/e_phi_max*sign(a_phi2);
P.kd_phi = 0.041; %(2*zeta_phi*wn_phi-a_phi1)/a_phi2;

% rlocus(tf([a_phi2],[1,(a_phi1+a_phi2*P.kd_phi),a_phi2*P.kp_phi,0]))

P.ki_phi = 2.69; %.03;

%% Course Hold Gains
Vg = Va;

%%%DESIGN PARAMETERS%%%
zeta_chi = .9;
WX = 37;
%%%END%%%%%%%%%%%%%%%%%

wn_chi = wn_phi/WX; % bandwidth separation

% P.kp_chi = 1.53; %2*zeta_chi*wn_chi*Vg/P.gravity;
% P.ki_chi = 0.6885; %wn_chi^2*Vg/P.gravity;

P.kp_chi = 1.2;
P.ki_chi = 0.041;

P.kp_phi = 0.88;
P.kd_phi = 0.16;
P.ki_phi = 0.304;

%% Sideslip Hold Gains
delta_r_max = 45*pi/180;

%%%DESIGN PARAMETERS%%%%%%
e_beta_max = 10*pi/180;
zeta_beta = .7;
%%%END%%%%%%%%%%%%%%%%%%%%

a_beta1 = -P.rho*Va*P.S_wing/(2*P.mass)*P.C_Y_beta;
a_beta2 = P.rho*Va*P.S_wing/(2*P.mass)*P.C_Y_delta_r;

P.kp_beta = delta_r_max/e_beta_max*sign(a_beta2);
P.ki_beta = 1/a_beta2*((a_beta1 + a_beta2*P.kp_beta)/2/zeta_beta)^2;

%% Pitch Attitude Hold Gains
delta_e_max = 30*pi/180;

%%%DESIGN PARAMETERS%%%%%%
e_theta_max = 20*pi/180;
zeta_theta = 0.9;
%%%END%%%%%%%%%%%%%%%%%%%%

a_theta1 = -P.rho*Va^2*P.c*P.S_wing/(2*P.Jy)*P.C_m_q*P.c/(2*Va);
a_theta2 = -P.rho*Va^2*P.c*P.S_wing*P.C_m_alpha/(2*P.Jy);
a_theta3 = P.rho*Va^2*P.c*P.S_wing/(2*P.Jy)*P.C_m_delta_e;

P.kp_theta = delta_e_max/e_theta_max*sign(a_theta3);
wn_theta = sqrt(a_theta2 + delta_e_max/e_theta_max*abs(a_theta3));
P.kd_theta = (2*zeta_theta*wn_theta - a_theta1)/a_theta3;

RL_kptheta = tf(a_theta3,[1,(a_theta1+P.kd_theta*a_theta3),a_theta2]);
RL_kdtheta = tf([a_theta3,0],[1,a_theta1,a_theta2 + P.kp_theta*a_theta3]);

P.K_theta_DC = (P.kp_theta*a_theta3)/(a_theta2 + P.kp_theta*a_theta3);

%% Altitude Hold Gains

%%%DESIGN PARAMETERS%%%%%%
zeta_h = 0.7;
W_h = 30;
%%%END%%%%%%%%%%%%%%%%%%%%

wn_h = 1/W_h*wn_theta;

P.ki_h = wn_h^2/(P.K_theta_DC*Va);
P.kp_h = 2*zeta_h*wn_h/(P.K_theta_DC*Va);

%% Airspeed Hold with Pitch

%%%DESIGN PARAMETERS%%%%%%
W_V2 = 15;
zeta_V2 = .85;
%%%END%%%%%%%%%%%%%%%%%%%%

wn_V2 = wn_theta/W_V2;

delta_t_trim = u_trim(4);
delta_e_trim = u_trim(1);

Va_trim = sqrt(x_trim(4)^2 + x_trim(5)^2 + x_trim(6)^2);
alpha_trim = atan(x_trim(6)/x_trim(4));

a_V1 = P.rho*Va*P.S_wing/P.mass*(P.C_D_0 + P.C_D_alpha*alpha_trim + ...
    P.C_D_delta_e*delta_e_trim) + P.rho*P.S_prop*P.C_prop*Va/P.mass;
a_V2 = P.rho*P.S_prop*P.C_prop*P.k_motor^2*delta_t_trim/P.mass;
a_V3 = P.gravity*cos(x_trim(8) - x_trim(9));

P.ki_V2 = -wn_V2^2/(P.K_theta_DC*P.gravity);
P.kp_V2 = (a_V1 - 2*zeta_V2*wn_V2)/(P.K_theta_DC*P.gravity);

%% Airspeed Hold with Throttle

%%%DESIGN PARAMETERS%%%%%%
wn_V = 4;
zeta_V = 0.95;
%%%END%%%%%%%%%%%%%%%%%%%%

P.ki_V = wn_V^2/a_V2;
P.kp_V = (2*zeta_V*wn_V - a_V1)/a_V2;
