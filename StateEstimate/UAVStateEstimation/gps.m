% gps.m
%   Compute the output of gps sensor
%

function y = gps(uu, P)

    % relabel the inputs
    Va      = uu(1);
%    alpha   = uu(2);
%    beta    = uu(3);
    wn      = uu(4);
    we      = uu(5);
%    wd      = uu(6);
    pn      = uu(7);
    pe      = uu(8);
    pd      = uu(9);
%    u       = uu(10);
%    v       = uu(11);
%    w       = uu(12);
%    phi     = uu(13);
%    theta   = uu(14);
    psi     = uu(15);
%    p       = uu(16);
%    q       = uu(17);
%    r       = uu(18);
    t       = uu(19);
    
 
    % construct North, East, and altitude GPS measurements
    % old measurements
    persistent v_n_d1;
    persistent v_e_d1;
    persistent v_d_d1;
        
    if t == 0
        v_n_d1 = 0;
        v_e_d1 = 0;
        v_d_d1 = 0;
    end        
    
    v_n = exp(-1/1100)*v_n_d1 + 0.4*randn(1);
    v_e = exp(-1/1100)*v_e_d1 + 0.4*randn(1);
    v_d = exp(-1/1100)*v_d_d1 + 0.7*randn(1);
    
    y_gps_n = pn + v_n;
    y_gps_e = pe + v_e; 
    y_gps_h = pd + v_d; 
    
    % update old error
    v_n_d1 = v_n;
    v_e_d1 = v_e;
    v_d_d1 = v_d;
    
    % construct groundspeed and course measurements
    Vn = Va*cos(psi) + wn;
    Ve = Va*sin(psi) + we;
    Vg = sqrt(Vn^2 + Ve^2);
    
    sigma_Vg = 2.1;
    sigma_chi = sigma_Vg/Vg;
    
    y_gps_Vg = Vg + sigma_Vg*randn(1);
    y_gps_course = atan2(Ve,Vn) + sigma_chi*randn(1);

    % construct total output
    y = [...
        y_gps_n;...
        y_gps_e;...
        y_gps_h;...
        y_gps_Vg;...
        y_gps_course;...
        ];
    
end



