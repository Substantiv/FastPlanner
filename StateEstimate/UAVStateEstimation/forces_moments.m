function out = forces_moments(x, delta, wind, P)

% Output :-
%       F - forces
%       M - moments
%       Va - airspeed
%       alpha - angle of attack
%       beta - sideslip
%       wind - wind vector in inertial frame

% Organize inputs
% states
pn = x(1);
pe = x(2);
pd = x(3);
u = x(4);
v = x(5);
w = x(6);
phi = x(7);
theta = x(8);
psi = x(9);
p = x(10);
q = x(11);
r = x(12);
% control surface positions
delta_e = delta(1);
delta_a = delta(2);
delta_r = delta(3);
delta_t = delta(4);
% steady wind
w_ns = wind(1);
w_es = wind(2);
w_ds = wind(3);
% gusts
u_wg = wind(4);
v_wg = wind(5);
w_wg = wind(6);
%

% compute wind data in NED
u_w = cos(theta)*cos(psi)*w_ns + cos(theta)*sin(psi)*w_es ...
    - sin(theta)*w_ds + u_wg;
v_w = (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi))*w_ns + ...
    (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi))*w_es + ...
    sin(phi)*cos(theta)*w_ds + v_wg;
w_w = (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*w_ns + ...
    (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w_es + ...
    cos(phi)*cos(theta)*w_ds + w_wg;

% compute airspeed in body frame
u_r = u - u_w;
v_r = v - v_w;
w_r = w - w_w;

% Airspeed magnitudes and angles
Va = sqrt(u_r^2 + v_r^2 + w_r^2);
alpha = atan(w_r/u_r);
beta = asin(v_r/sqrt(u_r^2 + v_r^2 + w_r^2));


if (Va == 0)
    Va = P.Va;
end
if ~isfinite(alpha)
    alpha = 0;
end
if ~isfinite(beta)
    beta = 0;
end

% transform wind back into inertial frame for output
w_n = cos(theta)*cos(psi)*u_w + (sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))*v_w...
    + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*w_w;
w_e = cos(theta)*sin(psi)*u_w + (sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi))*v_w...
    + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*w_w;
w_d = -sin(theta)*u_w + sin(phi)*cos(theta)*v_w + cos(phi)*cos(theta)*w_w;

% Coefficients
sigma = (1 + exp(-P.M*(alpha-P.alpha0)) + exp(P.M*(alpha+P.alpha0)))/...
    ((1 + exp(-P.M*(alpha-P.alpha0)))*(1 + exp(P.M*(alpha+P.alpha0))));
C_L_alpha = (1-sigma)*(P.C_L_0 + P.C_L_alpha*alpha) + sigma*(2*sign(alpha)*...
    sin(alpha)^2*cos(alpha));
C_D_alpha = P.C_D_p + (P.C_L_0 + P.C_L_alpha*alpha)^2/(pi*P.e*P.AR);

C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
C_X_q_alpha = -P.C_D_q*cos(alpha) + P.C_L_q*sin(alpha);
C_X_delta_e_alpha = -P.C_D_delta_e*cos(alpha) + P.C_L_delta_e*sin(alpha);
C_Z_alpha = -C_D_alpha*sin(alpha) - C_L_alpha*cos(alpha);
C_Z_q_alpha = -P.C_D_q*sin(alpha) - P.C_L_q*cos(alpha);
C_Z_delta_e_alpha = -P.C_D_delta_e*sin(alpha) - P.C_L_delta_e*cos(alpha);

  % compute external forces and torques on aircraft
  % Forces
fx = -P.mass*P.gravity*sin(theta) + 1/2*P.rho*Va^2*P.S_wing*(C_X_alpha + ...
    C_X_q_alpha*P.c*q/(2*Va) + C_X_delta_e_alpha*delta_e) + ...
    1/2*P.rho*P.S_prop*P.C_prop*((P.k_motor*delta_t)^2 - Va^2);
fy = P.mass*P.gravity*cos(theta)*sin(phi) + 1/2*P.rho*Va^2*P.S_wing*(P.C_Y_0 + ...
    P.C_Y_beta*beta + P.C_Y_p*P.b*p/(2*Va) + P.C_Y_r*P.b*r/(2*Va) + ...
    P.C_Y_delta_a*delta_a + P.C_Y_delta_r*delta_r);
fz = P.mass*P.gravity*cos(theta)*cos(phi) + 1/2*P.rho*Va^2*P.S_wing*(C_Z_alpha +...
    C_Z_q_alpha*P.c*q/(2*Va) + C_Z_delta_e_alpha*delta_e);

% Torques
ell = 1/2*P.rho*Va^2*P.S_wing*P.b*(P.C_ell_0 + P.C_ell_beta*beta + ...
    P.C_ell_p*P.b*p/(2*Va) + P.C_ell_r*P.b*r/(2*Va) + P.C_ell_delta_a*delta_a ...
    + P.C_ell_delta_r*delta_r) - P.k_T_P*(P.k_Omega*delta_t)^2;
m = 1/2*P.rho*Va^2*P.S_wing*P.c*(P.C_m_0 + P.C_m_alpha*alpha + ...
    P.C_m_q*P.c*q/(2*Va) + P.C_m_delta_e*delta_e);
n = 1/2*P.rho*Va^2*P.S_wing*P.b*(P.C_n_0 + P.C_n_beta*beta + ...
    P.C_n_p*P.b*p/(2*Va) + P.C_n_r*P.b*r/(2*Va) + P.C_n_delta_a*delta_a ...
    + P.C_n_delta_r*delta_r);

% output vector
out = [fx; fy; fz; ell; m; n; Va; alpha; beta; w_n; w_e; w_d];

end