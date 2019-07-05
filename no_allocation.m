%this function can transform velocity vector in body-fixed frame 
%to which in earth-fixed frame

%the Input varity 't' is a 12*1 matrix:
%t(1) ~ t(6)  Force in 6 axis;
%t(6) ~ t(12) Eta_Eth 状态.
function N = no_allocation(t)

T = [t(1);t(2);t(3);t(4);t(5);t(6)];
phi = t(10); theta = t(11); psi = t(12);

J = [cos(psi)*cos(theta),-sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi),sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta),0,0,0;
     sin(psi)*cos(theta),cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi),-cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi),0,0,0;
      -sin(theta)       ,                cos(theta)*sin(phi)           ,                      cos(theta)*cos(phi)      ,0,0,0;
                  0     ,                               0              ,                               0               ,1,sin(phi)*tan(theta),cos(phi)*tan(theta);
                  0     ,                               0              ,                               0               ,0,   cos(phi)        ,-sin(phi);
                  0     ,                               0              ,                               0               ,0,sin(phi)/cos(theta),cos(phi)/cos(theta)];
            
av =  0.45; % 垂直方向推进器x坐标绝对值
bv =  0.3; % 垂直方向推进器y坐标绝对值
% cv = -0.1; % 垂直方向推进器z坐标
cv = 0;
ah =  0.5; % 水平方向推进器x坐标绝对值
bh =  0.2; % 水平方向推进器y坐标绝对值
% ch =  0.2; % 水平方向推进器z坐标
ch = 0;
alpha = pi/6; %水平方向推进器与x轴夹角30度
%% 垂直方向4个推进器
e1 = [0,0,1]';      e7 = [0,0,1]';      e4 = [0,0,1]';      e6 = [0,0,1]'; % 推力向量
d1 = [-av,-bv, cv]';d7 = [ av,-bv, cv]';d4 = [-av, bv, cv]';d6 = [ av, bv, cv]'; % 方向向量
%% 水平方向4个推进器
e5 = [-cos(alpha), sin(alpha), 0]';  e3 = [ cos(alpha), sin(alpha), 0]'; 
e8 = [-cos(alpha),-sin(alpha), 0]';  e2 = [ cos(alpha),-sin(alpha), 0]'; % 推力向量
d5 = [-ah,-bh, ch]';                d3 = [ ah,-bh, ch]';
d8 = [-ah, bh, ch]';                d2 = [ ah, bh, ch]'; % 方向向量

%% 力和力矩向量
V1 = [e1;cross(d1,e1)]; V2 = [e2;cross(d2,e2)]; % 力和力矩向量
V3 = [e3;cross(d3,e3)]; V4 = [e4;cross(d4,e4)];
V5 = [e5;cross(d5,e5)]; V6 = [e6;cross(d6,e6)];
V7 = [e7;cross(d7,e7)]; V8 = [e8;cross(d8,e8)];
 
%% 推力配置矩阵
K = [V1 V2 V3 V4 V5 V6 V7 V8]; %推力配置矩阵
P = eye(8); %推进器使能矩阵

TEMP = K*[50 50 -50 -50 0 0 0 0]';
KMAX = abs(TEMP(4)); disp('KMAX ='); %disp(KMAX);%KMAX
TEMP = K*[50 -50 50 -50 0 0 0 0]';
MMAX = abs(TEMP(5)); disp('MMAX =');% disp(MMAX);%MMAX
TEMP = K*[0 0 0 0 50 -50 -50 50]';
NMAX = abs(TEMP(6)); disp('NMAX =');% disp(NMAX);%NMAX
TEMP = K*[0 0 0 0 50 -50 50 -50]';
XMAX = abs(TEMP(1)); disp('NMAX ='); %disp(NMAX);%XMAX
TEMP = K*[0 0 0 0 50 50 -50 -50]';
YMAX = abs(TEMP(2)); disp('NMAX =');% disp(NMAX);%YMAX
TEMP = K*[50 50 50 50 0 0 0 0]';
ZMAX = abs(TEMP(3)); disp('NMAX ='); %disp(NMAX);%ZMAX

if(T(1)>XMAX)
    T(1) = XMAX;
end

if(T(2)>YMAX)
    T(2) = YMAX;
end

if(T(3)>ZMAX)
    T(3) = ZMAX;
end

if(T(4)>KMAX)
    T(4) = KMAX;
end

if(T(5)>MMAX)
    T(5) = MMAX;
end

if(T(6)>NMAX)
    T(6) = NMAX;
end

%% 8推进器推力矩阵

fmin = -50*ones(8,1); % 推进器推力下限
fmax =  50*ones(8,1); % 推进器推力上限

%% 未分配算法

f = K'/J*T;

N = sign(f).*sqrt(abs(f)/1.926e-5);
end
