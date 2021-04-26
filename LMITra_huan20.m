% clear all
% clc

Lhuan20=...
  [  2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1;
    -1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1;
    -1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2;
   ];

[netG,L]=genNetbyL(Lhuan20);
lamda=[eig(L)]';
%lamda'=[-0.0000    0.0979    0.0979    0.3820    0.3820    0.8244    0.8244    1.3820    1.3820    2.0000
%  2.0000    2.6180    2.6180    3.1756    3.1756    3.6180    3.6180    3.9021    3.9021    4.0000
A=[1 1;0 1];
B=[1 1;0 1];

B2=-lamda(1,2)*B;B4=-lamda(1,3)*B;B6=-lamda(1,4)*B;B8=-lamda(1,5)*B;B10=-lamda(1,6)*B;B12=-lamda(1,7)*B;B14=-lamda(1,8)*B;B16=-lamda(1,9)*B;B18=-lamda(1,10)*B;B20=-lamda(1,11)*B;
B22=-lamda(1,12)*B;B24=-lamda(1,13)*B;B26=-lamda(1,14)*B;B28=-lamda(1,15)*B;B30=-lamda(1,16)*B;B32=-lamda(1,17)*B;B34=-lamda(1,18)*B;B36=-lamda(1,19)*B;B38=-lamda(1,20)*B;
% Bnew(1,2)=-lamda(1,1)*B;
% for i=3:1:20
%     j=2*(i-2);
% Bnew(j,j)=-lamda(i).*B
% end

% lamda1=0;lamda2=0.0698;lamda3=0.1390;lamda4=0.2679;lamda5=0.2679;lamda6=0.5967;lamda7=0.9344;lamda8=1.0000;lamda9=1.0000;lamda10=1.0000;lamda11=1.0000;lamda12=1.7853;lamda13=2.2057;
% lamda14=2.6042;lamda15=3.3312;lamda16=3.7321;lamda17=3.7321;lamda18=4.3475;lamda19=4.8322;lamda20=5.1541;
% 
%  B2=-lamda2*B;B4=-lamda3*B;B6=-lamda4*B;B8=-lamda5*B;B10=-lamda6*B;B12=-lamda7*B;B14=-lamda8*B;B16=-lamda9*B;B18=-lamda10*B;B20=-lamda11*B;
% B22=-lamda12*B;B24=-lamda13*B;B26=-lamda14*B;B28=-lamda15*B;B30=-lamda16*B;B32=-lamda17*B;B34=-lamda18*B;B36=-lamda19*B;B38=-lamda20*B;

% B2=-0.5858*B; 
% B4=-2*B;
% B6=-3.4142*B;


R=10^-16*[1 0;0 1];%0416过程噪声
Q=10^-14*[1 1;1 2];%测量噪声

% Q=10^-12*[1 0;0 1];%改了！！！！！！！！0413
% R=1.6*10^-11*[1 1;1 2];
% Q=10^-6*[1 0;0 1];
% R=0*[1 1;1 2];
l1=[1;1;1;1];
l=eye(4)-1/4*l1*l1';

% U=[-0.6533 -0.500 -0.2706;-0.2706 0.500 0.6533;0.2706 0.500 -0.6533;0.6533 -0.5000 0.2706];%这是哪里算出来的�?

R1=A*Q*A';
R2=R+R1+Q;
K=[0.1 0;0 0.1];

taoe=10^-10;
k2=1;%定义的一个中间变�?
taok=10^-10;
sz=1000;
sz1=1000;
% % sz=100;
% % sz1=100;
G=0;G2=0;
J=zeros(1,sz1);




setlmis([]);
P = lmivar(1, [2 1]);
P2 = lmivar(1, [2 1]);
P3 = lmivar(1, [2 1]);
P4 = lmivar(1, [2 1]);P5 = lmivar(1, [2 1]);P6 = lmivar(1, [2 1]);
P7 = lmivar(1, [2 1]);P8 = lmivar(1, [2 1]);P9 = lmivar(1, [2 1]);P10 = lmivar(1, [2 1]);P11 = lmivar(1, [2 1]);P12 = lmivar(1, [2 1]);
 P13 = lmivar(1, [2 1]);P14 = lmivar(1, [2 1]);P15 = lmivar(1, [2 1]);P16 = lmivar(1, [2 1]);P17 = lmivar(1, [2 1]);P18 = lmivar(1, [2 1]);
P19 = lmivar(1, [2 1]);
X  = lmivar(2,[2 2]);
H  = lmivar(2,[2 2]);

lmiterm([-4 1 1 X],1,1);

lmiterm([-1 1 1 P],1,1);%p-h
lmiterm([-1 1 1 0],-R2);
lmiterm([-1 1 2 X],A,1);%AX+B2H
lmiterm([-1 1 2 H],B2,1);
lmiterm([-1 2 2 X],1,1,'s');%X+X'-P
lmiterm([-1 2 2 P],-1,1);

lmiterm([-2 1 1 P2],1,1);
lmiterm([-2 1 1 0],-R2);
lmiterm([-2 1 2 X],A,1);
lmiterm([-2 1 2 H],B4,1);
lmiterm([-2 2 2 X],1,1,'s');
lmiterm([-2 2 2 P2],-1,1);

lmiterm([-3 1 1 P3],1,1);
lmiterm([-3 1 1 0],-R2);
lmiterm([-3 1 2 X],A,1);
lmiterm([-3 1 2 H],B6,1);
lmiterm([-3 2 2 X],1,1,'s');
lmiterm([-3 2 2 P3],-1,1);

lmiterm([-4 1 1 P4],1,1);
lmiterm([-4 1 1 0],-R2);
lmiterm([-4 1 2 X],A,1);
lmiterm([-4 1 2 H],B8,1);
lmiterm([-4 2 2 X],1,1,'s');
lmiterm([-4 2 2 P4],-1,1);

lmiterm([-5 1 1 P5],1,1);
lmiterm([-5 1 1 0],-R2);
lmiterm([-5 1 2 X],A,1);
lmiterm([-5 1 2 H],B10,1);
lmiterm([-5 2 2 X],1,1,'s');
lmiterm([-5 2 2 P5],-1,1);

lmiterm([-6 1 1 P6],1,1);
lmiterm([-6 1 1 0],-R2);
lmiterm([-6 1 2 X],A,1);
lmiterm([-6 1 2 H],B12,1);
lmiterm([-6 2 2 X],1,1,'s');
lmiterm([-6 2 2 P6],-1,1);

lmiterm([-7 1 1 P7],1,1);
lmiterm([-7 1 1 0],-R2);
lmiterm([-7 1 2 X],A,1);
lmiterm([-7 1 2 H],B14,1);
lmiterm([-7 2 2 X],1,1,'s');
lmiterm([-7 2 2 P7],-1,1);

lmiterm([-8 1 1 P8],1,1);
lmiterm([-8 1 1 0],-R2);
lmiterm([-8 1 2 X],A,1);
lmiterm([-8 1 2 H],B16,1);
lmiterm([-8 2 2 X],1,1,'s');
lmiterm([-8 2 2 P8],-1,1);

lmiterm([-9 1 1 P9],1,1);
lmiterm([-9 1 1 0],-R2);
lmiterm([-9 1 2 X],A,1);
lmiterm([-9 1 2 H],B18,1);
lmiterm([-9 2 2 X],1,1,'s');
lmiterm([-9 2 2 P9],-1,1);

lmiterm([-10 1 1 P10],1,1);
lmiterm([-10 1 1 0],-R2);
lmiterm([-10 1 2 X],A,1);
lmiterm([-10 1 2 H],B20,1);
lmiterm([-10 2 2 X],1,1,'s');
lmiterm([-10 2 2 P10],-1,1);

lmiterm([-11 1 1 P11],1,1);
lmiterm([-11 1 1 0],-R2);
lmiterm([-11 1 2 X],A,1);
lmiterm([-11 1 2 H],B22,1);
lmiterm([-11 2 2 X],1,1,'s');
lmiterm([-11 2 2 P11],-1,1);

lmiterm([-12 1 1 P12],1,1);
lmiterm([-12 1 1 0],-R2);
lmiterm([-12 1 2 X],A,1);
lmiterm([-12 1 2 H],B24,1);
lmiterm([-12 2 2 X],1,1,'s');
lmiterm([-12 2 2 P12],-1,1);

lmiterm([-13 1 1 P13],1,1);
lmiterm([-13 1 1 0],-R2);
lmiterm([-13 1 2 X],A,1);
lmiterm([-13 1 2 H],B26,1);
lmiterm([-13 2 2 X],1,1,'s');
lmiterm([-13 2 2 P13],-1,1);

lmiterm([-14 1 1 P14],1,1);
lmiterm([-14 1 1 0],-R2);
lmiterm([-14 1 2 X],A,1);
lmiterm([-14 1 2 H],B28,1);
lmiterm([-14 2 2 X],1,1,'s');
lmiterm([-14 2 2 P14],-1,1);

lmiterm([-15 1 1 P15],1,1);
lmiterm([-15 1 1 0],-R2);
lmiterm([-15 1 2 X],A,1);
lmiterm([-15 1 2 H],B30,1);
lmiterm([-15 2 2 X],1,1,'s');
lmiterm([-15 2 2 P15],-1,1);


lmiterm([-16 1 1 P16],1,1);
lmiterm([-16 1 1 0],-R2);
lmiterm([-16 1 2 X],A,1);
lmiterm([-16 1 2 H],B32,1);
lmiterm([-16 2 2 X],1,1,'s');
lmiterm([-16 2 2 P16],-1,1);

lmiterm([-17 1 1 P17],1,1);
lmiterm([-17 1 1 0],-R2);
lmiterm([-17 1 2 X],A,1);
lmiterm([-17 1 2 H],B34,1);
lmiterm([-17 2 2 X],1,1,'s');
lmiterm([-17 2 2 P17],-1,1);

lmiterm([-18 1 1 P18],1,1);
lmiterm([-18 1 1 0],-R2);
lmiterm([-18 1 2 X],A,1);
lmiterm([-18 1 2 H],B36,1);
lmiterm([-18 2 2 X],1,1,'s');
lmiterm([-18 2 2 P18],-1,1);

lmiterm([-19 1 1 P19],1,1);
lmiterm([-1 1 1 0],-R2);
lmiterm([-1 1 2 X],A,1);
lmiterm([-1 1 2 H],B38,1);
lmiterm([-1 2 2 X],1,1,'s');
lmiterm([-1 2 2 P19],-1,1);


lmis = getlmis;
options=[1e-8 0 0 0 0];
% c=mat2dec(lmis,eye(2),eye(2),eye(2),zeros(2,2),zeros(2,2));
c=mat2dec(lmis,eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2));

[copt xopt]=mincx(lmis,c,options); 
P=dec2mat(lmis,xopt,P);
P2=dec2mat(lmis,xopt,P2);
P3=dec2mat(lmis,xopt,P3);
P4=dec2mat(lmis,xopt,P4);P5=dec2mat(lmis,xopt,P5);P6=dec2mat(lmis,xopt,P6);
P7=dec2mat(lmis,xopt,P7);P8=dec2mat(lmis,xopt,P8);P9=dec2mat(lmis,xopt,P9);P10=dec2mat(lmis,xopt,P10);P11=dec2mat(lmis,xopt,P11);
P12=dec2mat(lmis,xopt,P12);
P13=dec2mat(lmis,xopt,P13);P14=dec2mat(lmis,xopt,P14);P15=dec2mat(lmis,xopt,P15);P16=dec2mat(lmis,xopt,P16);P17=dec2mat(lmis,xopt,P17);P18=dec2mat(lmis,xopt,P18);
P19=dec2mat(lmis,xopt,P19)
J1=copt;
X=dec2mat(lmis,xopt,X);
H=dec2mat(lmis,xopt,H);
K=H*inv(X);
save K
