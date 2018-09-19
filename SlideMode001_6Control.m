function Result=SlideMode001_6Control(tt)
%参数：质量1，密度1，长度1，设定位移6，设定速度6，设定加速度6，实际位移6,实际速度6，浮力6,水动力6  
%本函数为ROV的六自由度滑膜控制器
%经仿真发现，三个平移方向的自由度控制效果良好，可以实现解耦(未通过理论证明)，但是旋转控制未能实现。
%通过调节C值可以改变控制效果，改变k、delta、eps未实验
%本算法中未按照完全将状态方程模型代入到控制率中的方法，而是有取舍的代入到控制率中

%以后的改进：可以适当调节参数C、k、delta、eps，获得更好的效果；可以改变控制率中的取舍量。
mass=tt(1);rou=tt(2);L=tt(3);
setPosition=tt(4:9);%设定位移
setSpeed=tt(10:15);%设定速度
setAcc=tt(16:21);%设定加速度
actualPosition=tt(22:27);%实际位移
actualSpeed=tt(28:33);%实际速度
Fwb=tt(34:39);%浮力
Fwater=tt(40:45);
J=tt(46:48);
C=tt(49:54);
C=reshape(C,6,1);
Jx=J(1);
Jy=J(2);
Jz=J(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%水动力系数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X  
  Xud=-1.5777e-03;Yvd=-3.0753e-02;Zwd=-3.0799e-02;
  Kpd=-7.5371e-06;Mqd=-1.5970e-03; Nrd=-1.6012e-03;
%%%%%%%%%%%%%%%%%%%%【完】水动力系数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_=0.5*rou*L^3*Xud;         aa=mass-A_;
B_=0.5*rou*L^3*Yvd;         bb=mass-B_;
G_=0.5*rou*L^3*Zwd;         cc=mass-G_;
H_=0.5*rou*L^5*Kpd;         dd=Jx-H_;
K_=0.5*rou*L^5*Mqd;         ee=Jy-K_;
M_=0.5*rou*L^5*Nrd;         ff=Jz-M_;

Fu=[aa;bb;cc;dd;ee;ff];

eps=0.001*[1;1;1;1;1;1];
k=10*[1;1;1;1;1;1];
delta=5*[1;1;1;1;1;1];
%C=[5;5;5;1;1;1];
%%%%%%%%%%%%%%水动力计算--缺项%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Caculate the Water Force
   
%%%%%%%%%%%%%%%%%%%【完】水动力计算--缺项%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=C.*(actualPosition-setPosition)+actualSpeed-setSpeed;
Result=Fu.*(-C.*(actualPosition-setPosition)-(Fwater+Fwb)./Fu+setAcc+(-eps.*s./(abs(s)+delta)-k.*s));
%Result=[Re(1);Re(2);Re(3);Re(4);Re(5);Re(6)];
end