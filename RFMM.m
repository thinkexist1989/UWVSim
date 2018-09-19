%this function can produce the Restoring Forces and Moments Matrix 
%(short for RFMM)of the Underwater Welding Vehicle

%the Input varity 't' is a 11*1 matrix:
%[B;x_b;y_b;z_b;W;x_g;y_g;z_g;phi;theta;psi]

function g = RFMM(t)
    %t = reshape(t,11,1);
    B = t(1);%Buoyancy of UWV
    x_b = t(2); y_b = t(3); z_b = t(4); %the Center of Buoyancy
    W = t(5); %Weightof UWV
    x_g = t(6); y_g = t(7); z_g = t(8); %the Center of Gravity
    phi = t(9); theta = t(10); psi = t(11); %Orientation of UWV
    
    g = [-(W-B)*sin(theta);
          (W-B)*cos(theta)*sin(phi);
          (W-B)*cos(theta)*cos(phi);
          (y_g*W-y_b*B)*cos(theta)*cos(phi)-(z_g*W-z_b*B)*cos(theta)*sin(phi);
         -(z_g*W-z_b*B)*sin(theta)-(x_g*W-x_b*B)*cos(theta)*cos(phi);
          (x_g*W-x_b*B)*cos(theta)*sin(phi)+(y_g*W-y_b*B)*sin(theta)];
end