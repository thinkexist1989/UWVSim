%this function can calculate the hydrodynamics of
%underwater welding vehicle
function tau_h = Hydrodynamics(t)
%the Input varity 't' is a 14*1 matrix:
%[u;v;w;p;q;r;du;dv;dw;dp;dq;dr;L;rou]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define Hydrodynamic Factor:These Factors are 无因次量
%t is a 1*14 matrix varity including Length of ROV ,Density of Water, u v w p q r
    %ud vd wd pd qd rd

%define Hydrodynamic Factor:These Factors are 无因次量
%X  
    Xud=-1.5777e-02;    Xpr=1.1004e-05; Xqq=-9.8510e-05;    Xrr=-9.4196e-05;
    Xvr=1.750018e-05;     Xuu=-4.7337e-03;     % Xwp=-3.0799e-02;  Xuq=-3.0799e-02;
%Y
    Yrd=-Xrr; Ypd=-Xpr;    Ypq=-Xqq;
    Yvd=-3.0753e-02;    Ywp=-3.0799e-02; Yr=1.097760e-02;
      Y0=0;   Yv=-4.4961e-02;     Yv_v_=-1.6687e-01;
    Yr_r_=1.2579e-02;   
   Yv_r_=-1.3816e-02;
%Z
     Zwd=Ywp;    Zvp=Yvd;   Zq=-1.7092e-02; Zw_q_=-1.0737e-02;
     Z0=1.7856e-03;  Zw=-3.0407e-02;    Zw_w_=-1.4567e-01;
     Zqd=-Ypq;   Zpp=Ypd;%Zq_q_=-9.7750e-03;
     Zrp=Yrd;    
%K
    Kpd=-7.5371e-04;    Krd=-1.5074e-08;    Kqr=-4.2208e-08;
    Kpq=-4.0914e-08;    Kvd=-1.7607e-07;    Kvq=-4.31e-07;
         Kwr=0;    Kwp=0;%Kwq=1.1004e-05;
%M
    Mq=-1.0205e-02; Mw_q_=1.9550e-06;   M0=3.5067e-07;  Mw=1.0279e-07;
    Mrp=1.5936e-06; Mq_q_=-1.1447e-07;  Mwd=-1.0349e-06;    Mvr=0;
    Mvp=0;    Mw_w_=-5.3684e-06;  Mqd=-1.5970e-02; %Mqd=-1.5970e-03
    Mpp=4.0914e-07;     Mrr=-4.0914e-07;
%N
    Nrd=-1.6012e-03;    Npd=-4.0914e-06;    Npq=Mqd-Kpd;
    Nqr=-Npd;     Nr_r_=-1.4432e-03;  Nvd=1.0601e-03;
    Nwp=0;    Nr=-1.1679e-03;    % N_v_r=-2.5157e-02;
    N0=0;   Nv=-9.3775e-03;  Nvq=0;  % Nv_v_=6.9082e-03;    
       
%Other no defined factors

%Other no defined factors
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   % t=reshape(t,1,14);
%get all varity from t
    u = t(1); v = t(2); w = t(3);
    p = t(4); q = t(5); r = t(6);
    ud = t(7); vd = t(8); wd = t(9);
    pd = t(10); qd = t(11); rd = t(12);
    L = t(13);    %Length of ROV
    rou = t(14);  %Density of Water
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Caculate the Water Force
    X=0.5*rou*L^4*(Xqq*q^2+Xrr*r*r+Xpr*p*r)+...
        0.5*rou*L^3*(Xud*ud+Xvr*v*r)+...
        0.5*rou*L^2*(Xuu*u*u + Xud*u);
    
    Y=0.5*rou*L^4*(Yrd*rd+Ypd*pd+Yr_r_*r*abs(r))+...
        0.5*rou*L^3*(Yvd*vd+Ywp*w*p)+...
        0.5*rou*L^3*u*(Yr*r)+...
        0.5*rou*L^3*(Yv_r_*sign(v)*abs(sqrt(v*v+w*w))*abs(r))+...
        0.5*rou*L^2*(Y0*u*u+Yv*u*v+Yv_v_*v*abs(sqrt(v*v+w*w))+Yvd*v);
    
    Z=0.5*rou*L^4*(Zqd*qd+Zpp*p*p+Zrp*r*p)+...
        0.5*rou*L^3*(Zwd*wd+Zvp*v*p)+...
        0.5*rou*L^3*(Zq*u*q+Zw_q_*sign(w)*abs(sqrt(v*v+w*w))*abs(q))+...
        0.5*rou*L^2*(Z0*u*u+Zw*u*w+Zw_w_*w*abs(sqrt(v*v+w*w))+Zwd*w);
        
    
    K=0.5*rou*L^5*(Kpd*pd+Krd*rd+Kqr*q*r+Kpq*p*q)+...
        0.5*rou*L^4*(Kvd*vd+Kvq*v*q+Kwp*w*p+Kwr*w*r+Kpd*p);

%    M=0.5*rou*L^5*(Mqd*qd)+0.5*rou*L^4*(Mqd*q);
    M=0.5*rou*L^5*(Mqd*qd+Mpp*p*p+Mrr*r*r+Mrp*r*p+Mq_q_*q*abs(q))+...
        0.5*rou*L^4*(Mwd*wd+Mvr*v*r+Mvp*v*p+Mqd*q)+...
        0.5*rou*L^3*(Mq*u*q+Mw_q_*abs(sqrt(v*v+w*w))*q+M0*u*u+Mw*u*w+...
                     Mw_w_*w*abs(sqrt(v*v+w*w)));
                 
    N=0.5*rou*L^5*(Nrd*rd+Npd*pd+Npq*p*q+Nqr*q*r+Nr_r_*r*abs(r))+...
        0.5*rou*L^4*(Nvd*vd+Nwp*w*p+Nvq*v*q+Nrd*r)+...
        0.5*rou*L^3*(Nr*u*r+N0*u*u+Nv*u*v);   
    
   tau_h=[X;Y;Z;K;M;N];
end
    
        
