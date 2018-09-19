
phi = 0; theta = 0; psi = pi/2;
J = [cos(psi)*cos(theta)  -sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi)  sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta);
          sin(psi)*cos(theta)  cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi)   -cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi);
          -sin(theta)          cos(theta)*sin(phi)                              cos(theta)*cos(phi)];

      
n = [1 0 0]*J