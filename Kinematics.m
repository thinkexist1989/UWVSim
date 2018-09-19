%this function can transform velocity vector in body-fixed frame 
%to which in earth-fixed frame

%the Input varity 't' is a 9*1 matrix:
%[u;v;w;p;q;r;x;y;z;phi;theta;psi] 
%t(1:6,1) is the linear and angular velocity vector in body-fixed frame;
%t(7:12,1) is the position and orientation vector in earth-fixed frame.

function dEta = Kinematics(t)

    u = t(1,1);      v = t(2,1);        w = t(3,1); % linear velocity 
    p = t(4,1);      q = t(5,1);        r = t(6,1); % angular velocity
  %  x = t(7,1);      y = t(8,1);        z = t(9,1); % position
    phi = t(7,1);   theta = t(8,1);   psi = t(9,1);%orientation
    
    J1 = [cos(psi)*cos(theta)  -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi)  sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
          sin(psi)*cos(theta)  cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi)   -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
          -sin(theta)          cos(theta)*sin(phi)                              cos(theta)*cos(phi)];
    dEta_1 = J1*[u;v;w];
    
    J2 = [1  sin(phi)*tan(theta)   cos(phi)*tan(theta);
          0  cos(phi)              -sin(phi);
          0  sin(phi)/cos(theta)   cos(phi)/cos(theta);];
    dEta_2 = J2*[p;q;r];
    
    
    dEta = [dEta_1;dEta_2];
end