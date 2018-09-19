function Result=SMC(tt)
global deta_d_pre eta_d_pre deta_r_pre;



step = 0.1; %仿真步长

m = tt(1); % mass
xG = tt(2); yG = tt(3); zG = tt(4);        % center of gravity
    
Ix = tt(5);    Ixy = -tt(6);  Ixz = -tt(7); %inertia tensor matrix
Iyx = -tt(8);  Iy  =  tt(9);  Iyz = -tt(10);
Izx = -tt(11); Izy = -tt(12); Iz  =  tt(13);

g = tt(14:19);

x = tt(19);   y = tt(20);   z = tt(21);%position
phi = tt(22);   theta = tt(23);   psi = tt(24);%orientation
actualPosition = [x;y;z;phi;theta;psi];

u =tt(25);v =tt(26);w =tt(27);
p =tt(28);q =tt(29);r =tt(30); %state variable  uu is same as u 
actualSpeed = [u;v;w;p;q;r];

setPosition = tt(31:36); %设定设定位置

setSpeed = tt(37:42); %设定速度

M = [m     0      0      0       m*zG      -m*yG;   
     0     m      0    -m*zG      0         m*xG;
     0     0      m     m*yG    -m*xG          0;
     0     -m*zG  m*yG   Ix      -Ixy       -Ixz;
     m*zG  0      -m*xG  -Iyx     Iy        -Iyz;
    -m*yG  m*xG   0      -Izx    -Ixy         Iz];
    
C = [       0                  0                  0                 m*(yG*q+zG*r)         -m*(xG*q-w)           -m*(xG*r+v);
            0                  0                  0                -m*(yG*p+w)             m*(zG*r+xG*p)        -m*(yG*r-u);
            0                  0                  0                -m*(zG*p-v)            -m*(zG*q+u)           m*(xG*p+yG*q);
    -m*(yG*q + zG*r)       m*(yG*p+w)          m*(zG*p-v)               0                -Iyz*q-Ixz*p+Iz*r     Iyz*r+Ixy*p-Iy*q;
       m*(xG*q-w)         -m*(zG*r+xG*p)       m*(zG*q+u)         Iyz*q+Ixz*p-Iz*r           0                -Ixz*r-Ixy*q+Ix*p;
       m*(xG*r+v)          m*(yG*r-u)        -m*(xG*p+yG*q)      -Iyz*r-Ixy*p+Iy*q       Ixz*r+Ixy*q-Ix*p            0         ];

J1 = [cos(psi)*cos(theta)  -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi)  sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
      sin(psi)*cos(theta)  cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi)   -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
      -sin(theta)          cos(theta)*sin(phi)                              cos(theta)*cos(phi)];
    
J2 = [1  sin(phi)*tan(theta)   cos(phi)*tan(theta);
      0  cos(phi)              -sin(phi);
      0  sin(phi)/cos(theta)   cos(phi)/cos(theta);];

J = blkdiag(J1,J2);

error = actualPosition-setPosition;
derror = actualSpeed-setSpeed;
deta_d = (eta_d - eta_d_pre)/step;

deta_r = deta_d - C.*error;
ddeta_r = (deta_r - deta_r_pre)/step;

dalpha_r = J\deta_r;
ddalpha_r = J\ddeta_r;


%滑模面定义
s=derror+C.*error;
k = 0.1;
e = abs(M*ddalpha_r+C*dalpha_r+g);
tau = M*ddalpha_r+C*dalpha_r+g-e*sgn(J\s)-k*J\s;

deta_d_pre = deta_d;
eta_d_pre = eta_d;
deta_r_pre = deta_r;

Result = tau;
end