function [x_trim,u_trim] = compute_trim(filename, Va, gamma, R)
% Va is the desired airspeed (m/s)
% gamma is the desired flight path angle (radians)
% R is the desired radius (m) - use (+) for right handed orbit, 
%                                   (-) for left handed orbit


% define states
x0 = [0;...     %1 -  Pn
    0;...       %2 -  Pe
    0;...       %3 -  h
    Va;...      %4 -  u
    0;...       %5 -  v
    0;...       %6 -  w
    0;...       %7 -  phi
    gamma;...   %8 -  theta
    0;...       %9 -  psi
    0;...       %10-  p
    0;...       %11-  q
    0;...       %12-  r
    ];
% define which states to hold equal to the initial conditions
ix = [];

% specify initial inputs
u0 = [0;...     %1-  e
    0;...       %2-  a
    0;...       %3-  r
    1;...       %4-  t
    ];
% specify which inputs to hold constant
iu = [];

% define constant outputs
y0 = [Va;...    %1-  Va
    0;...       %2-  alpha
    0;...       %3-  beta
    ];
% specify which outputs to hold constant
iy = [1,3];

% define constant deerivatives
dx0 = [...              
    0;...               %1 - Pndot
    0;...               %2 - Pedot
    -Va*sin(gamma);...  %3 - hdot
    0;...               %4 - udot
    0;...               %5 - vdot
    0;...               %6 - wdot
    0;...               %7 - phidot
    0;...               %8 - thetadot
    0;...               %9 - psidot
    0;...               %10- pdot
    0;...               %11- qdot
    0;...               %12- rdot
    ];
if R~= 0, dx0(9) = Va/R; end
% specify which derivatives to hold constant in trim algrithm
idx = [3;4;5;6;7;8;9;10;11;12];


% compute trim conditions
[x_trim,u_trim,y_trim,dx_trim] = trim(filename,x0,u0,y0,ix,iu,iy,dx0,idx);

% check to make sure that the linearization worked (should be small)
norm(dx_trim(3:end)-dx0(3:end))

