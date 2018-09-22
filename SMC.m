%滑模控制算法
function Result=SMC(tt)

global error_sum;
global deta_d_pre;
global eta_d_pre;
global deta_r_pre;
global deta_pre;


if isempty(deta_d_pre)
    deta_d_pre = zeros(6,1);
end

if isempty(eta_d_pre)
    eta_d_pre = zeros(6,1);
end

if isempty(deta_r_pre)
    deta_r_pre = zeros(6,1);
end

if isempty(error_sum)
    error_sum = zeros(6,1);
end

if isempty(deta_pre)
    deta_pre = zeros(6,1);
end

step = 1; %仿真步长
cc = [0.1;0.1;0.1;0.1;0.1;0.1]; %滑模面参数

m = tt(1); % mass 1
xG = tt(2); yG = tt(3); zG = tt(4);        % center of gravity 2

Ix = tt(5);    Ixy = -tt(6);  Ixz = -tt(7); %inertia tensor matrix 3
Iyx = -tt(8);  Iy  =  tt(9);  Iyz = -tt(10);
Izx = -tt(11); Izy = -tt(12); Iz  =  tt(13);

g = tt(14:19); %4
g = reshape(g,6,1);

x = tt(20);   y = tt(21);   z = tt(22);%position  5
phi = tt(23);   theta = tt(24);   psi = tt(25);%orientation  5
actualPosition = [x;y;z;phi;theta;psi];

u =tt(26);v =tt(27);w =tt(28); %6
p =tt(29);q =tt(30);r =tt(31); %state variable  uu is same as u 6
actualSpeed = [u;v;w;p;q;r];

setPosition = tt(32:37); %设定设定位置 7
setPosition = reshape(setPosition,6,1);
eta_d = setPosition;

setSpeed = tt(38:43); %设定速度 8
setSpeed = reshape(setSpeed,6,1);


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

deta_r = deta_d - 1.*error;
ddeta_r = (deta_r - deta_r_pre)/step;

dalpha_r = J\deta_r;
ddalpha_r = J\ddeta_r;


%滑模面定义
s=cc.*derror+error;
k = 0.05*ones(6,1);
%e = 10*ones(6,1);
e = 3*abs(M*ddalpha_r+C*dalpha_r+g);
%tau = M*ddalpha_r+C*dalpha_r+g-e.*sign(J\s)-k.*(J\s);
%tau =M*ddalpha_r+C*dalpha_r-g-e.*sign(J\s)-k.*(J\s);
%tau = M*(actualSpeed-deta_pre)/0.1+C*actualSpeed-g-e.*sign(J\s)-k.*(J\s);
tau = M*(J\(actualSpeed-deta_pre))/0.1+C*(J\actualSpeed)-g-J*(e.*sign(J\s)-k.*(J\s));
deta_d_pre = deta_d;
eta_d_pre = eta_d;
deta_r_pre = deta_r;

error_sum = error_sum + error;
deta_pre = actualSpeed;

%Result = J*tau;
Result = tau;
end