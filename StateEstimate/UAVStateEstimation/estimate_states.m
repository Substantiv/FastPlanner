% estimate_states
%   - estimate the MAV states using gyros, accels, pressure sensors, and GPS.
%
% Outputs are:
%   pnhat    - estimated North position, 
%   pehat    - estimated East position, 
%   hhat     - estimated altitude, 
%   Vahat    - estimated airspeed, 
%   alphahat - estimated angle of attack
%   betahat  - estimated sideslip angle
%   phihat   - estimated roll angle, 
%   thetahat - estimated pitch angel, 
%   chihat   - estimated course angel, 
%   phat     - estimated roll rate, 
%   qhat     - estimated pitch rate, 
%   rhat     - estimated yaw rate,
%   Vghat    - estimated ground speed, 
%   wnhat    - estimate of North wind, 
%   wehat    - estimate of East wind
%   psihat   - estimate of heading angle


function xhat = estimate_states(uu,P)

    % rename inputs
    y_gyro_x      = uu(1);          % 陀螺仪测量x方向的角速度
    y_gyro_y      = uu(2);          % 陀螺仪测量y方向的角速度
    y_gyro_z      = uu(3);          % 陀螺仪测量z方向的角速度
    y_accel_x     = uu(4);          % 加速度计测量x方向的加速度
    y_accel_y     = uu(5);          % 加速度计测量y方向的加速度
    y_accel_z     = uu(6);          % 加速度计测量z方向的加速度
    y_static_pres = uu(7);          % 空速管测量的静压
    y_diff_pres   = uu(8);          % 空速管测量的差压
    y_gps_n       = uu(9);          % GPS测量的北位置
    y_gps_e       = uu(10);         % GPS测量的东位置
    y_gps_h       = uu(11);         % GPS测量的高度
    y_gps_Vg      = uu(12);         % GPS测量的地面速度
    y_gps_course  = uu(13);         % GPS测量的航向角
    t             = uu(14);         % 时间
    
    % Altitude and Airspeed estimation ------------------------------------
    % LPF inverted sensor models
    
    % Initialize values
    persistent y_static_pres_d1;
    persistent y_diff_pres_d1;

    if t == 0
       y_static_pres_d1 = 0;
       y_diff_pres_d1   = 0;
    end
   
    % Estimate angular rates
    lpf_static_pres = LPF(y_static_pres_d1,y_static_pres,P.alpha_static_pres);
    lpf_diff_pres   = LPF(y_diff_pres_d1,y_diff_pres,P.alpha_diff_pres);

    hhat  = lpf_static_pres/(P.rho*P.gravity);      % 高度估计值
    Vahat = sqrt(2*lpf_diff_pres/P.rho);            % 空速估计值
    phat  = y_gyro_x;                               % roll角速度估计值
    qhat  = y_gyro_y;                               % pitch角速度估计值
    rhat  = y_gyro_z;                               % yaw角速度估计值

    % update old values
    y_static_pres_d1 = lpf_static_pres;
    y_diff_pres_d1   = lpf_diff_pres;

    %----------------------------------------------------------------------
    % EKF Attitude Estimation
    %----------------------------------------------------------------------
    % Initialize states and covariance matrices
    persistent xhat_att;        % 姿态状态量的估计值
    persistent Q_att;           % 过程噪声的协方差
    persistent P_att;           % 估计误差的协方差
    persistent R_att;           % 测量噪声的协方差
    
    if t == 0
        xhat_att = [0, 0]';
        P_att = diag([0.01; 0]);
        Q_att = diag([1e-11, 1e-20]);
        R_att = 0.0025^2*eye(3);
    end
    
    % 姿态估计的输入
    u_att = [phat, qhat, rhat, Vahat];
    y_att = [y_accel_x, y_accel_y, y_accel_z]';
    
    % Prediction Steps
    N = 10;
    for i = 1:N
        xhat_att = xhat_att + P.Ts/N*f_att(xhat_att,u_att);
        A_att = df_dx_att(xhat_att,u_att);
        P_att = P_att + P.Ts/N*(A_att*P_att + P_att*A_att' + Q_att);
    end
            
    % Measurement Correction
    if ~mod(t,P.Ts_gyros)
        C_att = dh_dx_att(xhat_att,u_att,P.gravity);
        % 卡尔曼增益    
        L_att = P_att*C_att'*inv(R_att + C_att*P_att*C_att');       
        % 后验协方差
        P_att = (eye(2) - L_att*C_att)*P_att;                       
        % 后验状态估计
        xhat_att = xhat_att + L_att*(y_att - h_att(xhat_att,u_att,P.gravity));      
    end
    
    %----------------------------------------------------------------------
    % EKF GPS Estimation
    %----------------------------------------------------------------------
    % Initialize states and covariance matrices
    persistent xhat_gps;        % GPS状态量的估计值
    persistent Q_gps;           % 过程噪声的协方差
    persistent P_gps;           % 估计误差的协方差
    persistent R_gps;           % 测量噪声的协方差
    
    if t == 0
        xhat_gps = [P.pn0; P.pe0; P.Va; P.psi0; 0; 0; P.psi0];
        P_gps = diag([.01, .01, .01, .01, .01, .01, .01]);
        Q_gps = diag([.001, .001, .01, 0.01, .01, .01, .01]);
        R_gps = diag([5^2, 5^2, 2^2, (.45)^2, 20^2, 20^2]);
    end
    
    u_gps = [Vahat; qhat; rhat; xhat_att(1); xhat_att(2)];
    y_gps = [y_gps_n;...
            y_gps_e;...
            y_gps_Vg;...
            y_gps_course;...
            Vahat*cos(xhat_gps(7)) + xhat_gps(5) - xhat_gps(3)*cos(xhat_gps(4));...
            Vahat*sin(xhat_gps(7)) + xhat_gps(6) - xhat_gps(3)*sin(xhat_gps(4));...
            ];
    
    % Prediction Steps
    N = 10;
    for i = 1:N
        xhat_gps = xhat_gps + P.Ts/N*f_gps(xhat_gps,u_gps,P.gravity);
        A_gps = df_dx_gps(xhat_gps,u_gps,P.gravity);
        P_gps = P_gps + P.Ts/N*(A_gps*P_gps + P_gps*A_gps' + Q_gps);
    end
    
    % Measurement Correction
    if ~mod(t,P.Ts_gps)
        C_gps = dh_dx_gps(xhat_gps,u_gps);
        L_gps = P_gps*C_gps'*inv(R_gps + C_gps*P_gps*C_gps');
        P_gps = (eye(7) - L_gps*C_gps)*P_gps;
        xhat_gps = xhat_gps + L_gps*(y_gps - h_gps(xhat_gps,u_gps));
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % States out
    %----------------------------------------------------------------------
    % extract from attitude and gps estimations assign into xhat
    phihat   = xhat_att(1);
    thetahat = xhat_att(2);
    
    pnhat  = xhat_gps(1);
    pehat  = xhat_gps(2);
    Vghat  = xhat_gps(3);
    chihat = xhat_gps(4);
    wnhat  = xhat_gps(5);
    wehat  = xhat_gps(6);
    psihat = xhat_gps(7);
    
    alphahat = 0;
    betahat  = 0;
    bxhat    = 0;
    byhat    = 0;
    bzhat    = 0;
    
    xhat = [...
            pnhat;...
            pehat;...
            hhat;...
            Vahat;...
            alphahat;...
            betahat;...
            phihat;...
            thetahat;...
            chihat;...
            phat;...
            qhat;...
            rhat;...
            Vghat;...
            wnhat;...
            wehat;...
            psihat;...
            bxhat;...
            byhat;...
            bzhat;...
            ];
end

%----------------------------------------------------------------------
% Supporting Functions ************************************************
%----------------------------------------------------------------------
% Low-pass filter function: pass in previous filtered value, new 
% measurement, and filter gain value, get out new filtered value; 
% old value is updated to new one.
function y = LPF(y_d1, u, alpha)

    y = alpha*y_d1 + (1-alpha)*u;

end

% ATTITUDE FUNCTIONS *************************
% Attitude dynamics
function xhat_dot = f_att(xhat,u)
    p = u(1);
    q = u(2);
    r = u(3);
%     Va = u(4);
    phi = xhat(1);
    theta = xhat(2);
    
    xhat_dot = [...
                p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);...
                q*cos(phi) - r*sin(phi);...
                ];
end

% Jacobian of attitude dynamics
function A = df_dx_att(xhat,u)
%     p = u(1);
    q = u(2);
    r = u(3);
%     Va = u(4);
    phi = xhat(1);
    theta = xhat(2);
    
    A = [...
        q*cos(phi)*tan(theta)-r*sin(phi)*tan(theta), (q*sin(phi) + ...
        r*cos(phi))/cos(theta)^2;...

        -q*sin(phi) - r*cos(phi), 0;...
        ];
end

% State output for attitude
function h = h_att(xhat,u,g)
    p = u(1);
    q = u(2);
    r = u(3);
    Va = u(4);
    phi = xhat(1);
    theta = xhat(2);
    
    h = [...
        q*Va*sin(theta) + g*sin(theta);...
        r*Va*cos(theta) - p*Va*sin(theta) - g*cos(theta)*sin(phi);...
        -q*Va*cos(theta) - g*cos(theta)*cos(phi);...
        ];
end

% Jacobian of h for attitude
function C = dh_dx_att(xhat,u,g)
    p = u(1);
    q = u(2);
    r = u(3);
    Va = u(4);
    phi = xhat(1);
    theta = xhat(2);
    
    C = [...
        0, q*Va*cos(theta) + g*cos(theta);...
        
        -g*cos(theta)*cos(phi), -r*Va*sin(theta) - p*Va*cos(theta) + ...
        g*sin(theta)*sin(phi);...
        
        g*cos(theta)*sin(phi), q*Va*sin(theta) + g*sin(theta)*cos(phi);...
        ];
end

% GPS FUNCTIONS ********************************
% Nonlinear function for gps states
function xhat_dot = f_gps(xhat,u,g)

    % Enumerate states and inputs
%     pn = xhat(1);
%     pe = xhat(2);
    Vg = xhat(3);
    chi = xhat(4);
    wn = xhat(5);
    we = xhat(6);
    psi = xhat(7);
    
    Va = u(1);
    q = u(2);
    r = u(3);
    phi = u(4);
    theta = u(5);
    
    % nonlinear dynamics
    psidot = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);
    
    xhat_dot = [...
        Vg*cos(chi);...
        
        Vg*sin(chi);...
        
        ((Va*cos(psi) + wn)*(-Va*psidot*sin(psi)) + (Va*sin(psi) + we)...
        *(Va*psidot*cos(psi)))/Vg;...
        
        g/Vg*tan(phi)*cos(chi-psi);...
        
        0;...
        
        0;...
        
        psidot;...
        ];

end

% Jacobian of gps dynamics
function A = df_dx_gps(xhat,u,g)

    % Enumerate states and inputs
%     pn = xhat(1);
%     pe = xhat(2);
    Vg = xhat(3);
    chi = xhat(4);
    wn = xhat(5);
    we = xhat(6);
    psi = xhat(7);
    
    Va = u(1);
    q = u(2);
    r = u(3);
    phi = u(4);
    theta = u(5);
    
    psidot = q*sin(phi)/cos(theta) + r*cos(phi)/cos(theta);
    Vgdot = ((Va*cos(psi) + wn)*(-Va*psidot*sin(psi)) + (Va*sin(psi) + we)...
        *(Va*psidot*cos(psi)))/Vg;
    dVgdot_dpsi = -psidot*Va*(wn*cos(psi) + we*sin(psi))/Vg;
    dchidot_dVg = -g/Vg^2*tan(phi)*cos(chi-psi);
    dchidot_dchi = -g/Vg*tan(phi)*sin(chi-psi);
    dchidot_dpsi = g/Vg*tan(phi)*sin(chi-psi);
    
    % jacobian
    A = [...
        0, 0, cos(chi), -Vg*sin(chi), 0, 0, 0;...
        
        0, 0, sin(chi), Vg*cos(chi), 0, 0, 0;...
        
        0, 0, -Vgdot/Vg, 0, -psidot*Va*sin(psi)/Vg, ...
        psidot*Va*cos(psi)/Vg, dVgdot_dpsi;...
        
        0, 0, dchidot_dVg, dchidot_dchi, 0, 0, dchidot_dpsi;...
        
        0, 0, 0, 0, 0, 0, 0;...
        
        0, 0, 0, 0, 0, 0, 0;...
        
        0, 0, 0, 0, 0, 0, 0;...
        ];

end

% Output functions for gps model
function h = h_gps(xhat,u,g)
    % Enumerate states and inputs
    pn = xhat(1);
    pe = xhat(2);
    Vg = xhat(3);
    chi = xhat(4);
    wn = xhat(5);
    we = xhat(6);
    psi = xhat(7);
    
    Va = u(1);
%     q = u(2);
%     r = u(3);
%     phi = u(4);
%     theta = u(5);
    
    h = [pn;...
        pe;...
        Vg;...
        chi;...
        Va*cos(psi) + wn - Vg*cos(chi);...
        Va*sin(psi) + we - Vg*sin(chi);...
        ];
end

% Jacobian of gps sensor model
function C = dh_dx_gps(xhat,u)
    % Enumerate states and inputs
%     pn = xhat(1);
%     pe = xhat(2);
    Vg = xhat(3);
    chi = xhat(4);
%     wn = xhat(5);
%     we = xhat(6);
    psi = xhat(7);
    
    Va = u(1);
%     q = u(2);
%     r = u(3);
%     phi = u(4);
%     theta = u(5);
    
    C = [...
        1, 0, 0, 0, 0, 0, 0;...
        0, 1, 0, 0, 0, 0, 0;...
        0, 0, 1, 0, 0, 0, 0;...
        0, 0, 0, 1, 0, 0, 0;...
        0, 0, -cos(chi), Vg*sin(chi), 1, 0, -Va*sin(psi);...
        0, 0, -sin(chi), -Vg*cos(chi), 0, 1, Va*cos(psi);...
        ];
    
end
