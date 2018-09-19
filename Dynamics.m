%this function is s-function that calculating the dynamic model of
%underwater welding vehicle
function [sys,x0,str,ts,simStateCompliance] = Dynamics(t,x,u,flag)
global du dv dw dp dq dr;  %#ok<NUSED>
switch flag
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
  case 1
    sys=mdlDerivatives(t,x,u);
  case 2
    sys=mdlUpdate(t,x,u);
  case 3
    sys=mdlOutputs(t,x,u);
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);
  case 9
    sys=mdlTerminate(t,x,u);
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes
global du dv dw dp dq dr; 
du = 0; dv = 0; dw = 0; dp = 0; dq = 0; dr = 0;

sizes = simsizes;

sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 19;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

% initialize the initial conditions
x0  = [0;0;0;0;0;0]; % u v w p q r
% str is always an empty matrix
str = [];
% initialize the array of sample times
ts  = [0 0];

simStateCompliance = 'UnknownSimState';
% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%M*dv + C*v  = tau  tau = tau_h +tau_r + tau_p
function sys=mdlDerivatives(t,x,u)
global du dv dw dp dq dr; 
m = u(1); % mass
xG = u(2); yG = u(3); zG = u(4);        % center of gravity
    
Ix = u(5);    Ixy = -u(6);  Ixz = -u(7); %inertia tensor matrix
Iyx = -u(8);  Iy  =  u(9);  Iyz = -u(10);
Izx = -u(11); Izy = -u(12); Iz  =  u(13);
    
X = u(14); Y = u(15); Z = u(16);        % propeller thrusters
K = u(17); M = u(18); N = u(19);
    
uu =x(1);v =x(2);w =x(3);
p =x(4); q =x(5);r =x(6); %state variable  uu is same as u 

%tau = tau_h + tau_r + tau_p   
%tau_h is hydrodynamic 
%tau_r is restoring force and moments 
%tau_p is propeller's thrust
tau = [X;Y;Z;K;M;N]; 
Nu = [uu;v;w;p;q;r]; % state variable matrix

% udot = X/m + v*r - w*q + xG*(q^2 + r^2) - yG*(p*q - rdot) - zG*(p*r + qdot);
% vdot = Y/m + w*p - uu*r + yG*(r^2 + p^2) - zG*(q*r - pdot) - xG*(q*p + rdot);
% wdot = Z/m + uu*q - v*p + zG*(p^2 + q^2) - xG*(r*p  - qdot) - yG*(r*q +pdot);
% pdot = (K-(Iz-Iy)*q*r+(rdot+p*q)*Ixz-(r^2-q^2)*Iyz-(p*r-qdot)*Ixy-m*(yG*(wdot-uu*q+v*q)-zG*(vdot-w*p+uu*r)))/Ix;
% qdot = (M-(Ix-Iz)*r*p+(pdot+q*r)*Ixy-(p^2-r^2)*Izx-(q*p-rdot)*Iyz-m*(zG*(udot-v*r+w*q)-xG*(wdot-uu*q+v*p)))/Iy;
% rdot = (N-(Iy-Ix)*p*q+(qdot+r*p)*Iyz-(q^2-p^2)*Ixy-(r*q-pdot)*Izx-m*(xG*(vdot-w*p+uu*r)-yG*(udot-v*r+w*q)))/Iz;

%M
M = [m     0      0      0       m*zG      -m*yG;   
     0     m      0    -m*zG      0         m*xG;
     0     0      m     m*yG    -m*xG          0;
     0     -m*zG  m*yG   Ix      -Ixy       -Ixz;
     m*zG  0      -m*xG  -Iyx     Iy        -Iyz;
    -m*yG  m*xG   0      -Izx    -Ixy         Iz];
%C 
C = [       0                  0                  0                 m*(yG*q+zG*r)         -m*(xG*q-w)           -m*(xG*r+v);
            0                  0                  0                -m*(yG*p+w)             m*(zG*r+xG*p)        -m*(yG*r-uu);
            0                  0                  0                -m*(zG*p-v)            -m*(zG*q+uu)           m*(xG*p+yG*q);
     -m*(yG*q + zG*r)       m*(yG*p+w)          m*(zG*p-v)               0                -Iyz*q-Ixz*p+Iz*r     Iyz*r+Ixy*p-Iy*q;
        m*(xG*q-w)         -m*(zG*r+xG*p)       m*(zG*q+uu)         Iyz*q+Ixz*p-Iz*r           0                -Ixz*r-Ixy*q+Ix*p;
        m*(xG*r+v)          m*(yG*r-uu)        -m*(xG*p+yG*q)      -Iyz*r-Ixy*p+Iy*q       Ixz*r+Ixy*q-Ix*p            0         ];

dNu = inv(M)*(tau - C*Nu); %#ok<MINV>
%dNu = M\(tau - C*Nu); 

sys = dNu;
du = dNu(1); dv = dNu(2); dw = dNu(3); 
dp = dNu(4); dq = dNu(5); dr = dNu(6);

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
global du dv dw dp dq dr; 
sys = [x(1);x(2);x(3);x(4);x(5);x(6);du;dv;dw;dp;dq;dr];
%sys = [x(1);x(2);x(3);x(4);x(5);x(6)];
% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
