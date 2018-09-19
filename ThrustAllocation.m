%this function can transform velocity vector in body-fixed frame 
%to which in earth-fixed frame

%the Input varity 't' is a 12*1 matrix:
%t(1) ~ t(6)  Force in 6 axis;
%t(6) ~ t(12) Eta_Eth ״̬.
function N = ThrustAllocation(t)

T = [t(1);t(2);t(3);t(4);t(5);t(6)];
phi = t(10); theta = t(11); psi = t(12);

J = [cos(psi)*cos(theta),-sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi),sin(psi)*sin(phi)+cos(psi)*cos(phi)*sin(theta),0,0,0;
     sin(psi)*cos(theta),cos(psi)*cos(phi)+sin(phi)*sin(theta)*sin(psi),-cos(psi)*sin(phi)+sin(theta)*sin(psi)*cos(phi),0,0,0;
      -sin(theta)       ,                cos(theta)*sin(phi)           ,                      cos(theta)*cos(phi)      ,0,0,0;
                  0     ,                               0              ,                               0               ,1,sin(phi)*tan(theta),cos(phi)*tan(theta);
                  0     ,                               0              ,                               0               ,0,   cos(phi)        ,-sin(phi);
                  0     ,                               0              ,                               0               ,0,sin(phi)/cos(theta),cos(phi)/cos(theta)]
              
K_star =[-0.03     ,  -0.0882353 , 0.25  ,  -1.47059,    0.5  ,        0;
         0.288675  ,      -0.5   ,    0  ,     0    ,      0  , 0.629367;
         0.288675  ,    0.5      ,    0  ,    0     ,     0   ,-0.629367;
         -0.03     , 0.0882353   , 0.25  ,  1.47059 ,    0.5  ,        0;
         -0.288675 ,     0.5     ,  0    ,     0    ,      0  , 0.629367;
         0.03      , 0.0882353   ,  0.25 , 1.47059  ,   -0.5  ,        0;
         0.03      ,-0.0882353   ,  0.25 , -1.47059 ,   -0.5  ,        0;
        -0.288675  ,    -0.5     ,    0  ,     0    ,     0   , -0.629367];


f = K_star/J*T;

N = sign(f).*sqrt(abs(f)/1.926e-5);
end
