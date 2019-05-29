function Result=PID(tt)

persistent error_sum;
persistent deta_d_pre;
persistent eta_d_pre;
persistent deta_r_pre;
persistent deta_pre;


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

error = setPosition - actualPosition;
derror = setSpeed - actualSpeed;

J1 = [cos(psi)*cos(theta)  -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi)  sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
    sin(psi)*cos(theta)  cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi)   -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
    -sin(theta)          cos(theta)*sin(phi)                              cos(theta)*cos(phi)];

J2 = [1  sin(phi)*tan(theta)   cos(phi)*tan(theta);
      0  cos(phi)              -sin(phi);
      0  sin(phi)/cos(theta)   cos(phi)/cos(theta);];

J = blkdiag(J1,J2);

%M
M = [m     0      0      0       m*zG      -m*yG;   
     0     m      0    -m*zG      0         m*xG;
     0     0      m     m*yG    -m*xG          0;
     0     -m*zG  m*yG   Ix      -Ixy       -Ixz;
     m*zG  0      -m*xG  -Iyx     Iy        -Iyz;
    -m*yG  m*xG   0      -Izx    -Ixy         Iz];
%C 
C = [       0                  0                  0                 m*(yG*q+zG*r)         -m*(xG*q-w)           -m*(xG*r+v);
            0                  0                  0                -m*(yG*p+w)             m*(zG*r+xG*p)        -m*(yG*r-u);
            0                  0                  0                -m*(zG*p-v)            -m*(zG*q+u)           m*(xG*p+yG*q);
     -m*(yG*q + zG*r)       m*(yG*p+w)          m*(zG*p-v)               0                -Iyz*q-Ixz*p+Iz*r     Iyz*r+Ixy*p-Iy*q;
        m*(xG*q-w)         -m*(zG*r+xG*p)       m*(zG*q+u)         Iyz*q+Ixz*p-Iz*r           0                -Ixz*r-Ixy*q+Ix*p;
        m*(xG*r+v)          m*(yG*r-u)        -m*(xG*p+yG*q)      -Iyz*r-Ixy*p+Iy*q       Ixz*r+Ixy*q-Ix*p            0         ];


%%
%PID
P = [1;1;7;15;.06;.3];%P = [1;1;1;15;.06;.3];%P = [1;1;1;15;.06;.3]; 20
D = [0;0;0;12;0;1];    %D = [0;0;0;8;0;1];%D = [0;0;0;3;0;1]; 30
I = [0;0;0;0;0;0];    %I = [0;0;0;0;0;0];%I = [0;0;0;0.02;0;0];
tau = P.*(error)+I.*error_sum+D.*(derror);
%%
%Model-based PID
% P = -10*ones(6,1);
% D = -2*ones(6,1);
% I = -0.3*ones(6,1);
% tau = J*M/J*(actualSpeed-deta_pre)/0.1+J*C/J*actualSpeed-g+P.*(error)+I.*error_sum+D.*(derror);
%%
error_sum = error_sum + error;
deta_pre = actualSpeed;
Result = tau;

end