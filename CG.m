%this function can calculate the center of Gravity of UVW 

%the Input varity 't' is a 19*1 matrix:
%[X;Y;Z;xa;xb;xc;ya;yb;yc;za;zb;zc;mx;my;mz;iy;jy;ky;mr]

function rG_G = CG(t)
      
    X = t(1,1);   Y = t(2,1);   Z = t(3,1);   %Order
    xa = t(4,1);  xb = t(5,1);  xc = t(6,1);  %length width height of X axis
    ya = t(7,1);  yb = t(8,1);  yc = t(9,1);  %length width height of Y axis
    za = t(10,1); zb = t(11,1); zc = t(12,1); %length width height of Z axis
    mx = t(13,1); my = t(14,1); mz = t(15,1); %mass of X,Y,Z axises
    iy = t(16,1); jy = t(17,1); ky = t(18,1); %the coordinate of the CG of Y axis
    mr = t(19,1);                             %mass of UVW'body    
    
    g = 9.8;     %acceleration of gravity
    ryG = [iy;jy;ky]; % the coordinate of the CG of Y axis
    rxG = [iy;Y+jy;ky+(xc+yc)/2]; % the coordinate of the CG of X axis
    rzG = [X+iy;Y+jy;ky+xc+(yc+zc)/2]; % the coordinate of the CG of Y axis
    
    rG = (mx*rxG +my*ryG + mz*rzG)/(mr+mx+my+mz);
    G = mr*g + mx*g + my*g + mz*g;
    
    rG_G = [rG;G];   
end