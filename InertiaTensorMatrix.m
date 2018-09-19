%this function can calculate the Inertia Tensor Matrix 
%of the Three-Axises Platform
%the Input varity 't' is a 18*1 matrix:
%[X;Y;Z;xa;xb;xc;ya;yb;yc;za;zb;zc;mx;my;mz;iy;jy;ky]

function I = InertiaTensorMatrix(t)
      
    X = t(1);   Y = t(2);   Z = t(3);   %Order
    xa = t(4);  xb = t(5);  xc = t(6);  %length width height of X axis
    ya = t(7);  yb = t(8);  yc = t(9);  %length width height of Y axis
    za = t(10); zb = t(11); zc = t(12); %length width height of Z axis
    mx = t(13); my = t(14); mz = t(15); %mass of X,Y,Z axises
    iy = t(16); jy = t(17); ky = t(18); %the coordinate of the CG of Y axis
    
    ryG = [iy;jy;ky]; % the coordinate of the CG of Y axis
    rxG = [iy;Y+jy-yb/2;ky+(xc+yc)/2]; % the coordinate of the CG of X axis
    rzG = [X+iy-xa/2;Y+jy-yb/2;ky+xc+(yc+zc)/2]; % the coordinate of the CG of Y axis
    
    Ixc = InertiaTensorAboutCenter(mx,xa,xb,xc); %Rotation about CG
    Iyc = InertiaTensorAboutCenter(my,ya,yb,yc); %Rotation about CG
    Izc = InertiaTensorAboutCenter(mz,za,zb,zc); %Rotation about CG
    
    Ix = Ixc-mx*(rxG*rxG'-rxG'*rxG*eye(3)); 
    Iy = Iyc-my*(ryG*ryG'-ryG'*ryG*eye(3)); 
    Iz = Izc-mz*(rzG*rzG'-rzG'*rzG*eye(3)); 
    
    I = Ix+Iy+Iz;   
end

function Ic = InertiaTensorAboutCenter(m,a,b,c)
v = [m*(b^2+c^2)/12; m*(a^2+c^2)/12; m*(a^2 +b^2)/12];
Ic = diag(v,0);
end