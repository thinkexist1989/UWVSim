%this function can calculate the center of Buoyancy of UVW 

%the Input varity 't' is a 19*1 matrix:
%[X;Y;Z;xa;xb;xc;ya;yb;yc;za;zb;zc;mx;my;mz;iy;jy;ky;mr]

function rB_B = CB(t)
      
    X = t(1,1);   Y = t(2,1);   Z = t(3,1);   %Order
    xa = t(4,1);  xb = t(5,1);  xc = t(6,1);  %length width height of X axis
    ya = t(7,1);  yb = t(8,1);  yc = t(9,1);  %length width height of Y axis
    za = t(10,1); zb = t(11,1); zc = t(12,1); %length width height of Z axis
    mx = t(13,1); my = t(14,1); mz = t(15,1); %mass of X,Y,Z axises
    iy = t(16,1); jy = t(17,1); ky = t(18,1); %the coordinate of the CG of Y axis
    mr = t(19,1);                             %mass of UVW'body
    
    ryB = [iy;jy;ky]; % the coordinate of the CG of Y axis
    rxB = [iy;Y+jy;ky+(xc+yc)/2]; % the coordinate of the CG of X axis
    rzB = [X+iy;Y+jy;ky+xc+(yc+zc)/2]; % the coordinate of the CG of Y axis
    
    Rho = 1000;  %water density
    g = 9.8;     %acceleration of gravity
    Vx = xa*xb*xc; Vy = xa*xb*xc; Vz = xa*xb*xc; %volume of each axis
    Bx = Rho*g*Vx; By = Rho*g*Vy; Bz = Rho*g*Vz; %buoyancy of each axis
    
    Bp = mr*g + mx*g + my*g + mz*g - Bx - By - Bz;
    
    rB = (Bx*rxB + By*ryB +Bz*rzB)/(Bp + Bx +By +Bz);
    
    B = mr*g + mx*g + my*g + mz*g;
    
    rB_B = [rB;B];
    
end