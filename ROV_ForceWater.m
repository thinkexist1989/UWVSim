function result = ForceWater(t)
A = 0;
if A == 0
    %t is a 1*14 matrix varity including Length of ROV ,Density of Water, u v w p q r
    %ud vd wd pd qd rd

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define Hydrodynamic Factor:These Factors are 无因次量
    %X  
        Xdu = -213.4e-03;  Xdq = -1.134e-2;     %Xdu = -0.9996;
        Xuu = -2.626e-2;  Xvv = 0;    Xww = 0;    Xqq = -2.127e-2;  Xrr = 2.23e-2;
        Xuw = 0; 
        Xvr = 9.676e-3;   Xwq = -9.243e-3;  Xpr = 1.1004e-05;

    %Y
        Ydv = -305.2e-03;  Ydr = -1.663e-4;   Ydp = 1.63e-3;   %Ydv = -0.9676   Ydr = -0.223
        Yv = 2.1545e-2;   Yp = 1.47e-2;     Yr = -1.5e-3;  %Yv = -2.1545
        Yv_v_= -2.02e-2; Yp_p_ = 0;      Yr_r_ = -6.8e-3;
        Yvw = 0;        Yv_q_ = 0;      Yv_r_ = -2.04e-2;
        Yw_p_ = 9.243e-2; Yw_r_ = 0;
        Ypq = 9.851e-4;   Yqr = 0;   %Ypq = 0.2127

    %Z
        Zdw = -488.9e-3;  Zdq = -9.851e-4;     % Zdw = -0.9243    Zdq = -0.2127;
        Zw = -7.847e-2;   Zq = -5.38e-2;   Zvv = 0;
        Zww = -9.0e-3;   Zpp = 0;    Zqq = -1.98e-2;  Zrr = 0;
        Z_v_w = 0;     Zvp =-305.2e-04;  Zvr = 0;    Zw_q_ = 6.509e-4;   %Zvp =-0.9676
        Zpr = -1.663e-3;   %Zpr = -0.223
    %K
        Kdp = -4.1e-3;   Kdv = Ydp;
        Kv = 3.27e-2;    Kp = -1.32e-2;    Kr = 0;
        Kv_v_ = 0;      Kp_p_ = 0;      Kr_r_ = 0;
        Kvw = 0;        Kvq = -4.357e-5;  Kwp = 1.98e-2;   Kwr = -Ydr-Zdq;   %Kvq = -0.4357; Kwr = 0.4357
        Kpq = 0;        Kqr = -4e-3;

    %M
        Mdq = -1.37e-3;   Mdw = Zdq;  Mdu = Xdq;
        Mw = 3.3e-3;    Mq = -2.812e-3;   %Mq = -2.812;  
        Mvv = 0;        Mww = 3.2545e-4;   Mpp = 0;    Mq_q_ = 2.21e-5; Mrr = 0;
        M_v_w = 0;      Mvp = -Ydr;    Mvr = -1.98e-2;      M_w_q = 3.852e-4;     Mpr = 5.6e-3; %Mvp = 0.223;
    %N
        Ndr = 9.7e-3;   Ndv = Ydr;
        Nv = -1.5e-3;   Nr = -1.3696e-2;   Nv_v_ = 1.01e-2;      Np_p_ = 0;     Nr_r_ = -5.1e-3;    %Nr = -1.3696;
        Nvw = 0;        Nvq = 1.98e-2;   N_v_r = -1.36e-2;     Nwp = Zdq;  %Nwp = -0.2127;
        Npq =9.6e-3;     Nqr = 0;   

    %Other no defined factors

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        t=reshape(t,1,14);
    %get all varity from t    
        L=t(1,1);%Length of ROV
        rou=t(1,2);%Density of Water
        u =t(1,3);
        v =t(1,4);
        w =t(1,5);
        p =t(1,6);
        q =t(1,7);
        r =t(1,8);
        du =t(1,9);
        dv =t(1,10);
        dw =t(1,11);
        dp =t(1,12);
        dq =t(1,13);
        dr =t(1,14);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Caculate the Water Force
       X=0.5*rou*L^4*(Xqq*q*q+Xrr*r*r+Xpr*p*r)+...
            0.5*rou*L^3*(Xdu*du+Xvr*v*r+Xwq*w*q)+...
            0.5*rou*L^2*(Xuu*u*u+Xvv*v*v+Xww*w*w+Xuw*u*w);

        Y=0.5*rou*L^4*(Ydr*dr+Ydp*dp+Yp_p_*p*abs(p)+Yr_r_*r*abs(r)+Ypq*p*q+Yqr*q*r)+...
            0.5*rou*L^3*(Ydv*dv+Yw_p_*w*abs(p)+Yw_r_*w*abs(r))+...
            0.5*rou*L^3*u*(Yp*p+Yr*r)+...
            0.5*rou*L^3*(Yv_r_*sign(v)*abs(sqrt(v*v+w*w))*abs(r)+Yv_q_*sign(v)*abs(sqrt(v*v+w*w))*abs(q))+...
            0.5*rou*L^2*u*(Yv*v)+...
            0.5*rou*L^2*(Yv_v_*v*abs(sqrt(v*v+w*w))+Yvw*v*w);

        Z=0.5*rou*L^4*(Zdq*dq+Zpp*p*p+Zqq*q*q+Zrr*r*r+Zpr*p*r)+...
            0.5*rou*L^3*u*(Zq*q)+...
            0.5*rou*L^3*(Zdw*dw+Zvr*v*r+Zvp*v*p)+...
            0.5*rou*L^3*(Zw_q_*sign(w)*abs(sqrt(v*v+w*w))*abs(q))+...
            0.5*rou*L^2*u*(Zw*w)+...
            0.5*rou*L^2*(Zvv*v*v+Zww*abs(w*sqrt(v*v+w*w)));            %Z_v_w未加入


        K=0.5*rou*L^5*(Kdp*dp+Kqr*q*r+Kpq*p*q+Kp_p_*p*abs(p)+Kr_r_*r*abs(r))+...
            0.5*rou*L^4*u*(Kp*p+Kr*r)+...
            0.5*rou*L^4*(Kdv*dv+Kvq*v*q+Kwp*w*p+Kwr*w*r)+...
            0.5*rou*L^3*u*(Kv*v)+...
            0.5*rou*L^3*(Kv_v_*v*abs(sqrt(v*v+w*w)))+...
            0.5*rou*L^3*(Kvw*v*w);


        M=0.5*rou*L^5*(Mdq*dq+Mpp*p*p+Mrr*r*r+Mpr*p*r+Mq_q_*q*abs(q))+...
            0.5*rou*L^4*u*(Mq*q)+...
            0.5*rou*L^4*(Mdu*du+Mdw*dw+Mvr*v*r+Mvp*v*p)+...
            0.5*rou*L^3*u*(Mw*w)+...
            0.5*rou*L^3*(Mvv*v*v+M_w_q*abs(sqrt(v*v+w*w))*q+Mww*abs(w*sqrt(v*v+w*w))); %M_v_w未加入

        N=0.5*rou*L^5*(Ndr*dr+Npq*p*q+Nqr*q*r+Nr_r_*r*abs(r)+Np_p_*p*abs(p))+...
            0.5*rou*L^4*u*(Nr*r)+...
            0.5*rou*L^4*(Ndv*dv+Nwp*w*p+Nvq*v*q+N_v_r*abs(sqrt(v*v+w*w))*r)+...        
            0.5*rou*L^3*u*(Nv*v)+...
            0.5*rou*L^3*(Nvw*v*w+Nv_v_*v*abs(sqrt(v*v+w*w)));

elseif A == 1
    %t is a 1*14 matrix varity including Length of ROV ,Density of Water, u v w p q r
    %ud vd wd pd qd rd

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define Hydrodynamic Factor:These Factors are 无因次量
    %X  
        Xud=-1.5777e-03;    Xpr=1.1004e-05; Xqq=-9.8510e-04;    Xrr=-9.4196e-04;
        Xvr=1.750018E-02;     Xuu=-4.7337e-03;     % Xwp=-3.0799e-02;  Xuq=-3.0799e-02;
    %Y
        Yrd=-Xrr; Ypd=-Xpr;    Ypq=-Xqq;
        Yvd=-3.0753e-02;    Ywp=3.0799e-02; Yr=1.097760e-02;
          Y0=0;   Yv=-4.4961e-02;     Yv_v_=-1.6687e-01;
        Yr_r_=1.2579e-02;   
       Yv_r_=-1.3816e-02;
    %Z
         Zwd=-Ywp;    Zvp=Yvd;   Zq=-1.7092e-02; Zw_q_=-1.0737e-02;
         Z0=1.7856e-03;  Zw=-3.0407e-02;    Zw_w_=-1.4567e-01;
          Zqd=-Ypq;   Zpp=Ypd;%Zq_q_=-9.7750e-03;
         Zrp=Yrd;    
    %K
        Kpd=-7.5371e-06;    Krd=-1.5074e-06;    Kqr=-4.2208e-06;
        Kpq=-4.0914e-06;    Kvd=-1.7607e-05;    Kvq=-4.31e-05;
             Kwr=-Yrd-Zqd;    Kwp=-Ypd;%Kwq=1.1004e-05;
    %M
        Mq=-1.0205e-02; Mw_q_=1.9550e-02;   M0=3.5067e-04;  Mw=1.0279e-02;
        Mrp=1.5936e-03; Mq_q_=-1.1447e-03;  Mwd=-1.0349e-03;    Mvr=Ypd;
        Mvp=-Yrd;    Mw_w_=-5.3684e-03;  Mqd=-1.5970e-03; 
        Mpp=4.0914e-06;     Mrr=-4.0914e-06;
    %N
        Nrd=-1.6012e-03;    Npd=-4.0914e-06;    Npq=Mqd-Kpd;
        Nqr=-Npd;     Nr_r_=-1.4432e-03;  Nvd=1.0601e-03;
        Nwp=Zqd;    Nr=-1.1679e-02;    % N_v_r=-2.5157e-02;
        N0=0;   Nv=-9.3775e-03;  Nvq=-Ypd;  % Nv_v_=6.9082e-03;    

    %Other no defined factors

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        t=reshape(t,1,14);
    %get all varity from t
        L=t(1,1);%Length of ROV
        rou=t(1,2);%Density of Water
        u =t(1,3);
        v =t(1,4);
        w =t(1,5);
        p =t(1,6);
        q =t(1,7);
        r =t(1,8);
        ud =t(1,9);
        vd =t(1,10);
        wd =t(1,11);
        pd =t(1,12);
        qd =t(1,13);
        rd =t(1,14);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Caculate the Water Force
        X=0.5*rou*L^4*(Xqq*q^2+Xrr*r*r+Xpr*p*r)+...
            0.5*rou*L^3*(Xud*ud+Xvr*v*r)+...
            0.5*rou*L^2*(Xuu*u*u);

        Y=0.5*rou*L^4*(Yrd*rd+Ypd*pd+Yr_r_*r*abs(r))+...
            0.5*rou*L^3*(Yvd*vd+Ywp*w*p)+...
            0.5*rou*L^3*u*(Yr*r)+...
            0.5*rou*L^3*(Yv_r_*v/abs(v)*abs(sqrt(v*v+w*w))*abs(r))+...
            0.5*rou*L^2*(Y0*u*u+Yv*u*v+Yv_v_*v*abs(sqrt(v*v+w*w)));

        Z=0.5*rou*L^4*(Zqd*qd+Zpp*p*p+Zrp*r*p)+...
            0.5*rou*L^3*(Zwd*wd+Zvp*v*p)+...
            0.5*rou*L^3*(Zq*u*q+Zw_q_*w/abs(w)*abs(sqrt(v*v+w*w))*abs(q))+...
            0.5*rou*L^2*(Z0*u*u+Zw*u*w+Zw_w_*w*abs(sqrt(v*v+w*w)));


        K=0.5*rou*L^5*(Kpd*pd+Krd*rd+Kqr*q*r+Kpq*p*q)+...
            0.5*rou*L^4*(Kvd*vd+Kvq*v*q+Kwp*w*p+Kwr*w*r);


        M=0.5*rou*L^5*(Mqd*qd+Mpp*p*p+Mrr*r*r+Mrp*r*p+Mq_q_*q*abs(q))+...
            0.5*rou*L^4*(Mwd*wd+Mvr*v*r+Mvp*v*p)+...
            0.5*rou*L^3*(Mq*u*q+Mw_q_*abs(sqrt(v*v+w*w))*q+M0*u*u+Mw*u*w+...
                         Mw_w_*w*abs(sqrt(v*v+w*w)));

        N=0.5*rou*L^5*(Nrd*rd+Npd*pd+Npq*p*q+Nqr*q*r+Nr_r_*r*abs(r))+...
            0.5*rou*L^4*(Nvd*vd+Nwp*w*p+Nvq*v*q)+...
            0.5*rou*L^3*(Nr*u*r+N0*u*u+Nv*u*v);
end

result=[X;Y;Z;K;M;N];
end
    
        
