function y = autopilot(uu,P)
%
% autopilot for mavsim
% 

    % process inputs
    NN = 0;
%    pn       = uu(1+NN);  % inertial North position
%    pe       = uu(2+NN);  % inertial East position
    h        = uu(3+NN);  % altitude
    Va       = uu(4+NN);  % airspeed
    alpha    = uu(5+NN);  % angle of attack
    beta_air     = uu(6+NN);  % side slip angle
    phi      = uu(7+NN);  % roll angle
    theta    = uu(8+NN);  % pitch angle
    chi      = uu(9+NN);  % course angle
    p        = uu(10+NN); % body frame roll rate
    q        = uu(11+NN); % body frame pitch rate
    r        = uu(12+NN); % body frame yaw rate
%    Vg       = uu(13+NN); % ground speed
%    wn       = uu(14+NN); % wind North
%    we       = uu(15+NN); % wind East
%    psi      = uu(16+NN); % heading
%    bx       = uu(17+NN); % x-gyro bias
%    by       = uu(18+NN); % y-gyro bias
%    bz       = uu(19+NN); % z-gyro bias
    NN = NN+19;
    Va_c     = uu(1+NN);  % commanded airspeed (m/s)
    h_c      = uu(2+NN);  % commanded altitude (m)
    chi_c    = uu(3+NN);  % commanded course (rad)
    NN = NN+3;
    t        = uu(1+NN);   % time
    
    autopilot_version = 2;
        % autopilot_version == 1 <- used for tuning
        % autopilot_version == 2 <- standard autopilot defined in book
        % autopilot_version == 3 <- Total Energy Control for longitudinal AP
    switch autopilot_version
        case 1,
           [delta, x_command] = autopilot_tuning(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P);
        case 2,
           [delta, x_command] = autopilot_uavbook(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P);
        case 3,
               [delta, x_command] = autopilot_TECS(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P);
    end
    y = [delta; x_command];
end
    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autopilot versions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopilot_tuning
%   - used to tune each loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta, x_command] = autopilot_tuning(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P)

    mode = 2;
    switch mode
        case 1 % tune the roll loop
            phi_c = chi_c; % interpret chi_c to autopilot as course command
            
            flag = 0;
            
            if t == 0 
                flag = 1;
            end
            
            delta_a = roll_hold(phi_c,phi,p,flag,P);
            
            delta_r = 0; % no rudder
            % use trim values for elevator and throttle while tuning the lateral autopilot
            delta_e = P.u_trim(1);
            delta_t = P.u_trim(4);
            theta_c = 0;
        case 2, % tune the course loop
            if t==0,
                flag = 1;
            else
                flag = 0;
            end
            
            phi_c   = course_hold(chi_c, chi, r, flag, P);
            delta_a = roll_hold(phi_c, phi, p, flag, P);
            delta_r = 0; % no rudder
            % use trim values for elevator and throttle while tuning the lateral autopilot
            delta_e = P.u_trim(1);
            delta_t = P.u_trim(4);
            theta_c = 0;
        case 3, % tune the throttle to airspeed loop and pitch loop simultaneously
            theta_c = 20*pi/180; % + h_c;
            chi_c = 0;
            if t==0,
                
                phi_c   = course_hold(chi_c, chi, r, 1, P);
                delta_a = roll_hold(phi_c, phi, p, 1, P);

                delta_t = airspeed_throttle_hold(Va_c, Va, 1, P);
           else
                phi_c   = course_hold(chi_c, chi, r, 0, P);
                delta_a = roll_hold(phi_c, phi, p, 0, P);
                
                delta_t = airspeed_throttle_hold(Va_c, Va, 0, P);
            end
            
            delta_e = pitch_hold(theta_c, theta, q, P);
            
            delta_r = 0; % no rudder
            % use trim values for elevator and throttle while tuning the lateral autopilot
        case 4, % tune the pitch to airspeed loop 
            chi_c = 0;
            delta_t = P.u_trim(4);
            if t==0,
                phi_c   = course_hold(chi_c, chi, r, 1, P);
                theta_c = airspeed_pitch_hold(Va_c, Va, 1, P);
                delta_a = roll_hold(phi_c, phi, p, 1, P);
           else
                phi_c   = course_hold(chi_c, chi, r, 0, P);
                theta_c = airspeed_pitch_hold(Va_c, Va, 0, P);
                delta_a = roll_hold(phi_c, phi, p, 0, P);
            end
            
            delta_e = pitch_hold(theta_c, theta, q, P);
            delta_r = 0; % no rudder
            % use trim values for elevator and throttle while tuning the lateral autopilot
        case 5, % tune the pitch to altitude loop 
            chi_c = 0;
            if t==0,
                phi_c   = course_hold(chi_c, chi, r, 1, P);
                delta_a = roll_hold(phi_c, phi, p, 1, P);
                theta_c = altitude_hold(h_c, h, 1, P);
                delta_t = airspeed_throttle_hold(Va_c, Va, 1, P);
           else
                phi_c   = course_hold(chi_c, chi, r, 0, P);
                delta_a = roll_hold(phi_c, phi, p, 0, P);
                theta_c = altitude_hold(h_c, h, 0, P);
                delta_t = airspeed_throttle_hold(Va_c, Va, 0, P);
            end
            
            delta_e = pitch_hold(theta_c, theta, q, P);
            delta_r = 0; % no rudder
            % use trim values for elevator and throttle while tuning the lateral autopilot
      end
    %----------------------------------------------------------
    % create outputs
    
    % control outputs
    delta = [delta_e; delta_a; delta_r; delta_t];
    % commanded (desired) states
    x_command = [...
        0;...                    % pn
        0;...                    % pe
        h_c;...                  % h
        Va_c;...                 % Va
        0;...                    % alpha
        0;...                    % beta
        phi_c;...                % phi
        %theta_c*P.K_theta_DC;... % theta
        theta_c;
        chi_c;...                % chi
        0;...                    % p
        0;...                    % q
        0;...                    % r
        ];
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopilot_uavbook
%   - autopilot defined in the uavbook
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta, x_command] = autopilot_uavbook(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P)

    %----------------------------------------------------------
    % lateral autopilot
    if t==0

        delta_r = 0; %sideslip_hold(beta_air, 1, P);
        
        phi_c   = course_hold(chi_c, chi, r, 1, P);
        delta_a = roll_hold(phi_c, phi, p, 1, P);

    else
        phi_c   = course_hold(chi_c, chi, r, 0, P);
        delta_a = roll_hold(phi_c, phi, p, 0, P);
        
        delta_r = 0; %sideslip_hold(beta_air, 0, P);
    end
  
    
    %----------------------------------------------------------
    % longitudinal autopilot
    
    % define persistent variable for state of altitude state machine
    persistent altitude_state;
    
    % initialize persistent variable
    if t==0
                
        % initialize controllers
        airspeed_pitch_hold(Va_c, Va, 1, P);
        airspeed_throttle_hold(Va_c, Va, 1, P);
        altitude_hold(h_c, h, 1, P);
        
    end
    
    % start state machine for longitudinal control    
    if h<=P.altitude_take_off_zone    
            altitude_state = 1; % take-off zone
        elseif h<=h_c-P.altitude_hold_zone
            altitude_state = 2; % climb zone
        elseif h>=h_c+P.altitude_hold_zone
            altitude_state = 3; % descend zone
        else
            altitude_state = 4; % hold zone
    end
    
    % implement state machine
    switch altitude_state
        case 1  % in take-off zone
            delta_t = 1;
            theta_c = P.theta_take_off;
        case 2  % climb zone
             delta_t = 1;
             theta_c = airspeed_pitch_hold(Va_c, Va, 0, P);
        case 3 % descend zone
            delta_t = 0;
            theta_c = airspeed_pitch_hold(Va_c, Va, 0, P);
        case 4 % altitude hold zone
            delta_t = airspeed_throttle_hold(Va_c, Va, 0, P);
            theta_c = altitude_hold(h_c, h, q, P);
    end

    delta_e = pitch_hold(theta_c, theta, q, P);
    
    %----------------------------------------------------------
    % create outputs
    
    % control outputs
    delta = [delta_e; delta_a; delta_r; delta_t];
    % commanded (desired) states
    x_command = [...
        0;...                    % pn
        0;...                    % pe
        h_c;...                  % h
        Va_c;...                 % Va
        0;...                    % alpha
        0;...                    % beta
        phi_c;...                % phi
        %theta_c*P.K_theta_DC;... % theta
        theta_c;
        chi_c;...                % chi
        0;...                    % p
        0;...                    % q
        0;...                    % r
        ];
            
    y = [delta; x_command];
 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autopilot_TECS
%   - longitudinal autopilot based on total energy control systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta, x_command] = autopilot_TECS(Va_c,h_c,chi_c,Va,h,chi,phi,theta,p,q,r,t,P)

    %----------------------------------------------------------
    % lateral autopilot
    if t==0,
        % assume no rudder, therefore set delta_r=0
        delta_r = 0;%coordinated_turn_hold(beta, 1, P);
        phi_c   = course_hold(chi_c, chi, r, 1, P);

    else
        phi_c   = course_hold(chi_c, chi, r, 0, P);
        delta_r = 0;%coordinated_turn_hold(beta, 0, P);
    end
    delta_a = roll_hold(phi_c, phi, p, P);     
  
    
    %----------------------------------------------------------
    % longitudinal autopilot based on total energy control
    
    
    delta_e = 0;
    delta_t = 0;
 
    
    %----------------------------------------------------------
    % create outputs
    
    % control outputs
    delta = [delta_e; delta_a; delta_r; delta_t];
    % commanded (desired) states
    x_command = [...
        0;...                    % pn
        0;...                    % pe
        h_c;...                  % h
        Va_c;...                 % Va
        0;...                    % alpha
        0;...                    % beta
        phi_c;...                % phi
        %theta_c*P.K_theta_DC;... % theta
        theta_c;
        chi_c;...                % chi
        0;...                    % p
        0;...                    % q
        0;...                    % r
        ];
            
    y = [delta; x_command];
 
end
   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autopilot functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROLL *********************************************************
function delta_a = roll_hold(phi_c, phi, p, flag, P)
    persistent phi_integrator;
    persistent phi_error_d1;
    
    if flag == 1
        phi_integrator = 0;
        phi_error_d1 = 0;
    end
    
    error = phi_c - phi;
    
    % discrete integrator
    phi_integrator = phi_integrator + P.Ts/2*(error + phi_error_d1);
    
    % compute output command and saturate
    delta_a = sat(P.kp_phi*error + P.ki_phi*phi_integrator - P.kd_phi*p...
        ,P.delta_a_up,P.delta_a_down);
    
    phi_error_d1 = error; % store old error value
    
    % integrator anti-windup
    if P.ki_phi ~= 0
        u_unsat = P.kp_phi*error + P.ki_phi*phi_integrator + P.kd_phi*p;
        phi_integrator = phi_integrator + P.Ts/P.ki_phi*(delta_a - u_unsat);
    end
end

% COURSE HOLD***********************************************
function phi_c = course_hold(chi_c, chi, r, flag, P)
    persistent chi_integrator;
    persistent chi_error_d1;
    
    if flag == 1
        chi_integrator = 0;
        chi_error_d1 = 0;
    end
    
    error = chi_c - chi;
    
    % discrete integrator
    chi_integrator = chi_integrator + P.Ts/2*(error + chi_error_d1);
    
    % compute output command and saturate
    phi_c = P.kp_chi*error + P.ki_chi*chi_integrator;
    
    chi_error_d1 = error; % store old error value
    
%     % integrator anti-windup
    if P.ki_chi ~= 0
        u_unsat = P.kp_chi*error + P.ki_chi*chi_integrator;
        chi_integrator = chi_integrator + P.Ts/P.ki_chi*(phi_c - u_unsat);
    end
end

% SLIDESLIP HOLD************************************************
function delta_r = sideslip_hold(beta_air,flag,P)
    persistent beta_integrator;
    persistent beta_error_d1;
    
    if flag == 1
        beta_integrator = 0;
        beta_error_d1 = 0;
    end
    
    % discrete integrator
    beta_integrator = beta_integrator + P.Ts/2*(beta_air + beta_error_d1);
    
    % compute output command and saturate
    delta_r = sat(-P.kp_beta*beta_air - P.ki_beta*beta_integrator, P.delta_r_up,...
        P.delta_r_down);
    
    beta_error_d1 = beta_air; % store old error value
    
    % integrator anti-windup
    if P.ki_beta ~= 0
        u_unsat = P.kp_beta*beta_air + P.ki_beta*beta_integrator;
        beta_integrator = beta_integrator + P.Ts/P.ki_beta*(delta_r - u_unsat);
    end
end

% PITCH ALTITUDE HOLD************************************
function delta_e = pitch_hold(theta_c, theta, q, P)

    delta_e = sat(P.kp_theta*(theta_c - theta) - P.kd_theta*q,...
        P.delta_e_up,P.delta_e_down);

end

% ALTITUDE HOLD USING PITCH **********************************
function theta_c = altitude_hold(h_c, h, flag, P)
    persistent h_integrator;
    persistent h_error_d1;
    
    if flag == 1
        h_integrator = 0;
        h_error_d1 = 0;
    end
    
    error = h_c - h;
    
    % discrete integrator
    h_integrator = h_integrator + P.Ts/2*(error + h_error_d1);
    
    % compute output command
    theta_c = P.kp_h*error + P.ki_h*h_integrator;
    
    h_error_d1 = error; % store old error value
    
    % integrator anti-windup
    if P.ki_h ~= 0
        u_unsat = P.kp_h*error + P.ki_h*h_integrator;
        h_integrator = h_integrator + P.Ts/P.ki_h*(theta_c - u_unsat);
    end
end

% AIRSPEED HOLD USING PITCH ******************************************
function theta_c = airspeed_pitch_hold(Va_c, Va, flag, P)
    persistent Va_integrator;
    persistent Va_error_d1;
    
    if flag == 1
        Va_integrator = 0;
        Va_error_d1 = 0;
    end
    
    error = Va_c - Va;
    
    % discrete integrator
    Va_integrator = Va_integrator + P.Ts/2*(error + Va_error_d1);
    
    % compute output command and saturate
    theta_c = P.kp_V2*error + P.ki_V2*Va_integrator;
    
    Va_error_d1 = error; % store old error value
    
    % integrator anti-windup
    if P.ki_V2 ~= 0
        u_unsat = P.kp_V2*error + P.ki_V2*Va_integrator;
        Va_integrator = Va_integrator + P.Ts/P.ki_V2*(theta_c - u_unsat);
    end
end

% AIRSPEED HOLD USING THROTTLE ********************************
function delta_t = airspeed_throttle_hold(Va_c, Va, flag, P)
    persistent Va_integrator;
    persistent Va_error_d1;
    
    if flag == 1
        Va_integrator = 0;
        Va_error_d1 = 0;
    end
    
    error = Va_c - Va;
    
    % discrete integrator
    Va_integrator = Va_integrator + P.Ts/2*(error + Va_error_d1);
    
    % compute output command and saturate
    delta_t = sat(P.u_trim(4) + P.kp_V*error + P.ki_V*Va_integrator,1,0);
    
    Va_error_d1 = error; % store old error value
    
    % integrator anti-windup
    if P.ki_V ~= 0
        u_unsat = P.kp_V*error + P.ki_V*Va_integrator;
        Va_integrator = Va_integrator + P.Ts/P.ki_V*(delta_t - u_unsat);
    end

end



  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sat
%   - saturation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = sat(in, up_limit, low_limit)
  if in > up_limit
      out = up_limit;
  elseif in < low_limit
      out = low_limit;
  else
      out = in;
  end
end
  
 