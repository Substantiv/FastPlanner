P.gravity = 9.81;
   
%%%% Using Aerosonde UAV parameters.
%%%% physical parameters of airframe
P.mass = 25;
P.Jx   = 0.8244;
P.Jy   = 1.135;
P.Jz   = 1.759;
P.Jxz  = 0.1204;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% aerodynamic coefficients
P.S_wing        = 0.55;
P.b             = 2.8956;
P.c             = 0.18994;
P.S_prop        = 0.2027;
P.rho           = 1.2682;
P.k_motor       = 80;
P.k_T_P         = 0;
P.k_Omega       = 0;
P.e             = 0.9;
P.AR            = P.b^2/P.S_wing;

P.C_L_0         = 0.28;
P.C_L_alpha     = 3.45;
% P.C_L_alpha = pi*P.AR/(1+sqrt(1+(P.AR/2)^2)); % approximation given
P.C_L_q         = 0.0;
P.C_L_delta_e   = -0.36;
P.C_D_0         = 0.03;
P.C_D_alpha     = 0.30;
P.C_D_p         = 0.0437;
P.C_D_q         = 0.0;
P.C_D_delta_e   = 0.0;
P.C_m_0         = -0.02338;
P.C_m_alpha     = -0.38;
P.C_m_q         = -3.6;
P.C_m_delta_e   = -0.5;
P.C_Y_0         = 0.0;
P.C_Y_beta      = -0.98;
P.C_Y_p         = 0.0;
P.C_Y_r         = 0.0;
P.C_Y_delta_a   = 0.0;
P.C_Y_delta_r   = -0.17;
P.C_ell_0       = 0.0;
P.C_ell_beta    = -0.12;
P.C_ell_p       = -0.26;
P.C_ell_r       = 0.14;
P.C_ell_delta_a = 0.08;
P.C_ell_delta_r = 0.105;
P.C_n_0         = 0.0;
P.C_n_beta      = 0.25;
P.C_n_p         = 0.022;
P.C_n_r         = -0.35;
P.C_n_delta_a   = 0.06;
P.C_n_delta_r   = -0.032;
P.C_prop        = 1.0;
P.M             = 50;
P.epsilon       = 0.1592;
P.alpha0        = 0.4712;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gamma values
P.Gamma = P.Jx*P.Jz - P.Jxz^2;
P.Gamma1 = P.Jxz*(P.Jx - P.Jy + P.Jz)/P.Gamma;
P.Gamma2 = (P.Jz*(P.Jz - P.Jy) + P.Jxz^2)/P.Gamma;
P.Gamma3 = P.Jz/P.Gamma;
P.Gamma4 = P.Jxz/P.Gamma;
P.Gamma5 = (P.Jz - P.Jx)/P.Jy;
P.Gamma6 = P.Jxz/P.Jy;
P.Gamma7 = ((P.Jx - P.Jy)*P.Jx + P.Jxz^2)/P.Gamma;
P.Gamma8 = P.Jx/P.Gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% More coefficients
P.C_p_0 = P.Gamma3*P.C_ell_0 + P.Gamma4*P.C_n_0;
P.C_p_beta = P.Gamma3*P.C_ell_beta + P.Gamma4*P.C_n_beta;
P.C_p_p = P.Gamma3*P.C_ell_p + P.Gamma4*P.C_n_p;
P.C_p_r = P.Gamma3*P.C_ell_r + P.Gamma4*P.C_n_r;
P.C_p_delta_a = P.Gamma3*P.C_ell_delta_a + P.Gamma4*P.C_n_delta_a;
P.C_p_delta_r = P.Gamma3*P.C_ell_delta_r + P.Gamma4*P.C_n_delta_r;
P.C_r_0 = P.Gamma4*P.C_ell_0 + P.Gamma8*P.C_n_0;
P.C_r_beta = P.Gamma4*P.C_ell_beta + P.Gamma8*P.C_n_beta;
P.C_r_p = P.Gamma4*P.C_ell_p + P.Gamma8*P.C_n_p;
P.C_r_r = P.Gamma4*P.C_ell_r + P.Gamma8*P.C_n_r;
P.C_r_delta_a = P.Gamma4*P.C_ell_delta_a + P.Gamma8*P.C_n_delta_a;
P.C_r_delta_r = P.Gamma4*P.C_ell_delta_r + P.Gamma8*P.C_n_delta_r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% compute trim condition using 'mavsim_trim.slx'
P.Va = 35; 			% initial airspeed 				% in m/s
gamma = 0*pi/180; 	% initial flight path angle 	% in radians
R = inf; 			% initial turn radius 			% in m
h0    = 0;  % initial altitude

% autopilot sample rate
P.Ts = 0.01;
P.Ts_gps = .1;
P.Ts_gyros = .01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% first cut ai initial conditions
P.pn0    =  -1000; 	% initial North position
P.pe0    =  0; 		% initial East position
P.pd0    =  -h0; 		% initial Down position (negative altitude)
P.u0     =  P.Va; 	% initial velocity along body x-axis
P.v0     =  0; 		% initial velocity along body y-axis
P.w0     =  0; 		% initial velocity along body z-axis
P.phi0   =  0; 		% initial roll angle
P.theta0 =  0; 		% initial pitch angle
P.psi0   =  0; 		% initial yaw angle
P.p0     =  0; 		% initial body frame roll rate
P.q0     =  0; 		% initial body frame pitch rate
P.r0     =  0; 		% initial body frame yaw rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% wind parameters
P.wind_n = 0;
P.wind_e = 0;
P.wind_d = 0;
P.L_u = 200;
P.L_v = 200;
P.L_w = 50;
P.sigma_u = 1.06; 
P.sigma_v = 1.06;
P.sigma_w = .7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% run trim command
[x_trim, u_trim] = compute_trim('mavsim_trim',P.Va,gamma,R);
P.u_trim = u_trim;
P.x_trim = x_trim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set initial conditions to trim conditions
%%%% initial conditions
P.pn0    = -1000;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -100;  % initial Down position (negative altitude)
P.u0     = x_trim(4);  % initial velocity along body x-axis
P.v0     = x_trim(5);  % initial velocity along body y-axis
P.w0     = x_trim(6);  % initial velocity along body z-axis
P.phi0   = x_trim(7);  % initial roll angle
P.theta0 = x_trim(8);  % initial pitch angle
P.psi0   = x_trim(9);  % initial yaw angle
P.p0     = x_trim(10);  % initial body frame roll rate
P.q0     = x_trim(11);  % initial body frame pitch rate
P.r0     = x_trim(12);  % initial body frame yaw rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% compute different transfer functions
[T_phi_delta_a,T_chi_phi,T_theta_delta_e,T_h_theta,T_h_Va,T_Va_delta_t,T_Va_theta,T_v_delta_r]...
    = compute_tf_model(x_trim,u_trim,P);

%%%% linearize the equations of motion around trim conditions
[A_lon, B_lon, A_lat, B_lat] = compute_ss_model('mavsim_trim',x_trim,u_trim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Autopilot 
% Set control surface limits
P.delta_e_up = 45*pi/180;
P.delta_e_down = -45*pi/180;

P.delta_a_up = 30*pi/180;
P.delta_a_down = -30*pi/180;

P.delta_r_up = 45*pi/180;
P.delta_r_down = -45*pi/180;

P.delta_t_up = 1;
P.delta_t_down = 0;

% Take-off Values
P.altitude_hold_zone = 80;
P.altitude_take_off_zone = 50;
P.theta_take_off = 15*pi/180;

% Autopilot gains
run('compute_gains.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sensor Parameters
% Sensor Bias
P.bias_gyro_x = 0;
P.bias_gyro_y = 0;
P.bias_gyro_z = 0;

%%%% LPF Gain Values  % LPF - Low Pass Filter
P.alpha_static_pres = .50;
P.alpha_diff_pres = .20;
P.alpha_gyro_x = .25;
P.alpha_gyro_y = .25;
P.alpha_gyro_z = .25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
