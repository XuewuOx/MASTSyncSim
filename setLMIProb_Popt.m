% R3=R2+A1*A*Q*A'+A*Q*A'*A1'

%A*的奇数行是A-lamdaBK,A*的偶数行是A-1的转置
A1=A-B2*K;A2=A1';
A3=A-B4*K;A4=A3';
A5=A-B6*K;A6=A5';
A7=A-B8*K;A8=A7';
A9=A-B10*K;A10=A9';
A11=A-B12*K;A12=A11';
A13=A-B14*K;A14=A13';
A15=A-B16*K;A16=A15';
A17=A-B18*K;A18=A17';
A19=A-B20*K;A20=A19';
A21=A-B22*K;A22=A21';
A23=A-B24*K;A24=A23';
A25=A-B26*K;A26=A25';
A27=A-B28*K;A28=A27';
A29=A-B30*K;A30=A29';
A31=A-B32*K;A32=A31';
A33=A-B34*K;A34=A33';
A35=A-B36*K;A36=A35';
A37=A-B38*K;A38=A37';
% A39=A-B40*K;A40=A39';
% A41=A-B42*K;A42=A41';
% A43=A-B*K;A44=A43';
% A45=A-B*K;A46=A45';
% A47=A-B*K;A48=A47';



% Initial a LMI system
setlmis([]);%type=1定义块对角的对称矩阵，struct(r,1)表示第r个块的大小，struct(r,2)表示第r个块的类型，<1--全矩阵，0--标量，-1--零阵》
%P is a symmetric matrix,has a block size of 2 and this block and this
%block is symmetric 定义变量：
P = lmivar(1, [2 1]);P2 = lmivar(1, [2 1]);P3 = lmivar(1, [2 1]);P4 = lmivar(1, [2 1]);P5 = lmivar(1, [2 1]);P6 = lmivar(1, [2 1]);
P7 = lmivar(1, [2 1]);P8 = lmivar(1, [2 1]);P9 = lmivar(1, [2 1]);P10 = lmivar(1, [2 1]);P11 = lmivar(1, [2 1]);P12 = lmivar(1, [2 1]);
 P13 = lmivar(1, [2 1]);P14 = lmivar(1, [2 1]);P15 = lmivar(1, [2 1]);P16 = lmivar(1, [2 1]);P17 = lmivar(1, [2 1]);P18 = lmivar(1, [2 1]);
P19 = lmivar(1, [2 1]);% P20 = lmivar(1, [2 1]);P21 = lmivar(1, [2 1]);P22 = lmivar(1, [2 1]);P23 = lmivar(1, [2 1]);P24 = lmivar(1, [2 1]);
%描述变量的每一项 只需要描述上三角或者下三角元素，否则会描述成另一个LMI
% lmiterm(termID,A,B,flag)
% +p表示位于第p个线性矩阵不等式的左边，-p代表了这个项位于第p个不等式的右边。左边通常指较小的那边
% 对于外部变量，termID(2:3)取值为[0,0];对于左面和右边的内部变量位置为[i,j]
% 对于外部变量，termID(4)：取值为0;对于A*X*B,取值为X；对于A*X'*B,取值为-X
%flag选为s，A*X+X'A'可表示为lmiterm([1 1 1 -X],1,A')
lmiterm([-2 1 1 P], 1, 1);lmiterm([-5 1 1 P2], 1, 1);lmiterm([-6 1 1 P3], 1, 1);lmiterm([-2 1 1 P4], 1, 1);lmiterm([-5 1 1 P5], 1, 1);lmiterm([-6 1 1 P6], 1, 1);
lmiterm([-2 1 1 P7], 1, 1);lmiterm([-5 1 1 P8], 1, 1);lmiterm([-6 1 1 P9], 1, 1);lmiterm([-2 1 1 P10], 1, 1);lmiterm([-5 1 1 P11], 1, 1);
 lmiterm([-6 1 1 P12], 1, 1);
 lmiterm([-2 1 1 P13], 1, 1);lmiterm([-5 1 1 P14], 1, 1);lmiterm([-6 1 1 P15], 1, 1);lmiterm([-2 1 1 P16], 1, 1);lmiterm([-5 1 1 P17], 1, 1);lmiterm([-6 1 1 P18], 1, 1);
 lmiterm([-2 1 1 P19], 1, 1);%lmiterm([-5 1 1 P20], 1, 1);lmiterm([-6 1 1 P21], 1, 1);lmiterm([-2 1 1 P22], 1, 1);lmiterm([-5 1 1 P23], 1, 1);lmiterm([-6 1 1 P24], 1, 1);

% pos in (1, 1)
lmiterm([-1 1 1 P],1,1);lmiterm([-1 1 1 0],-R2);
% pos (1, 2)
lmiterm([-1 1 2 P], A1, 1);
% pos(2, 2)
lmiterm([-1 2 2 P], 1, 1);

% pos in (1, 1)
lmiterm([-2 1 1 P2],1,1);lmiterm([-2 1 1 0],-R2);
% pos (1, 2)
lmiterm([-2 1 2 P2], A2, 1);
% pos(2, 2)
lmiterm([-2 2 2 P2], 1, 1);

% pos in (1, 1)
lmiterm([-3 1 1 P3],1,1);lmiterm([-3 1 1 0],-R2);
% pos (1, 2)
lmiterm([-3 1 2 P3], A3, 1);
% pos(2, 2)
lmiterm([-3 2 2 P3], 1, 1);

% pos in (1, 1)
lmiterm([-4 1 1 P4],1,1);lmiterm([-4 1 1 0],-R2);
% pos (1, 2)
lmiterm([-4 1 2 P4], A4, 1);
% pos(2, 2)
lmiterm([-4 2 2 P4], 1, 1);

% pos in (1, 1)
lmiterm([-5 1 1 P5],1,1);lmiterm([-5 1 1 0],-R2);
% pos (1, 2)
lmiterm([-5 1 2 P5], A5, 1);
% pos(2, 2)
lmiterm([-5 2 2 P5], 1, 1);

% pos in (1, 1)
lmiterm([-6 1 1 P6],1,1);lmiterm([-6 1 1 0],-R2);
% pos (1, 2)
lmiterm([-6 1 2 P6], A6, 1);
% pos(2, 2)
lmiterm([-6 2 2 P6], 1, 1);

% pos in (1, 1)
lmiterm([-7 1 1 P7],1,1);lmiterm([-7 1 1 0],-R2);
% pos (1, 2)
lmiterm([-7 1 2 P7], A7, 1);
% pos(2, 2)
lmiterm([-7 2 2 P7], 1, 1);

% pos in (1, 1)
lmiterm([-8 1 1 P8],1,1);lmiterm([-8 1 1 0],-R2);
% pos (1, 2)
lmiterm([-8 1 2 P8], A8, 1);
% pos(2, 2)
lmiterm([-8 2 2 P8], 1, 1);

% pos in (1, 1)
lmiterm([-9 1 1 P9],1,1);lmiterm([-3 1 1 0],-R2);
% pos (1, 2)
lmiterm([-9 1 2 P9], A9, 1);
% pos(2, 2)
lmiterm([-9 2 2 P9], 1, 1);

% pos in (1, 1)
lmiterm([-10 1 1 P10],1,1);lmiterm([-10 1 1 0],-R2);
% pos (1, 2)
lmiterm([-10 1 2 P10], A10, 1);
% pos(2, 2)
lmiterm([-10 2 2 P10], 1, 1);

% pos in (1, 1)
lmiterm([-11 1 1 P11],1,1);lmiterm([-11 1 1 0],-R2);
% pos (1, 2)
lmiterm([-11 1 2 P11], A11, 1);
% pos(2, 2)
lmiterm([-11 2 2 P11], 1, 1);

% pos in (1, 1)
lmiterm([-12 1 1 P12],1,1);lmiterm([-12 1 1 0],-R2);
% pos (1, 2)
lmiterm([-12 1 2 P12], A12, 1);
% pos(2, 2)
lmiterm([-12 2 2 P12], 1, 1);

% pos in (1, 1)
lmiterm([-13 1 1 P13],1,1);lmiterm([-13 1 1 0],-R2);
% pos (1, 2)
lmiterm([-13 1 2 P13], A13, 1);
% pos(2, 2)
lmiterm([-13 2 2 P13], 1, 1);

% pos in (1, 1)
lmiterm([-14 1 1 P14],1,1);lmiterm([-14 1 1 0],-R2);
% pos (1, 2)
lmiterm([-14 1 2 P14], A14, 1);
% pos(2, 2)
lmiterm([-14 2 2 P14], 1, 1);

% pos in (1, 1)
lmiterm([-15 1 1 P15],1,1);lmiterm([-15 1 1 0],-R2);
% pos (1, 2)
lmiterm([-15 1 2 P15], A15, 1);
% pos(2, 2)
lmiterm([-15 2 2 P15], 1, 1);

% pos in (1, 1)
lmiterm([-16 1 1 P16],1,1);lmiterm([-16 1 1 0],-R2);
% pos (1, 2)
lmiterm([-16 1 2 P16], A16, 1);
% pos(2, 2)
lmiterm([-16 2 2 P16], 1, 1);

% pos in (1, 1)
lmiterm([-17 1 1 P17],1,1);lmiterm([-17 1 1 0],-R2);
% pos (1, 2)
lmiterm([-17 1 2 P17], A17, 1);
% pos(2, 2)
lmiterm([-17 2 2 P17], 1, 1);

% pos in (1, 1)
lmiterm([-18 1 1 P18],1,1);lmiterm([-18 1 1 0],-R2);
% pos (1, 2)
lmiterm([-18 1 2 P18], A18, 1);
% pos(2, 2)
lmiterm([-18 2 2 P18], 1, 1);

% pos in (1, 1)
lmiterm([-19 1 1 P19],1,1);lmiterm([-19 1 1 0],-R2);
% pos (1, 2)
lmiterm([-19 1 2 P19], A19, 1);
% pos(2, 2)
lmiterm([-19 2 2 P19], 1, 1);
% 
% % pos in (1, 1)
% lmiterm([-20 1 1 P20],1,1);lmiterm([-20 1 1 0],-R2);
% % pos (1, 2)
% lmiterm([-20 1 2 P20], A20, 1);
% % pos(2, 2)
% lmiterm([-20 2 2 P20], 1, 1);
% 
% % pos in (1, 1)
% lmiterm([-21 1 1 P21],1,1);lmiterm([-21 1 1 0],-R2);
% % pos (1, 2)
% lmiterm([-21 1 2 P21], A21, 1);
% % pos(2, 2)
% lmiterm([-21 2 2 P21], 1, 1);
% 
% % pos in (1, 1)
% lmiterm([-22 1 1 P22],1,1);lmiterm([-22 1 1 0],-R2);
% % pos (1, 2)
% lmiterm([-22 1 2 P22], A22, 1);
% % pos(2, 2)
% lmiterm([-22 2 2 P22], 1, 1);
% 
% % pos in (1, 1)
% lmiterm([-23 1 1 P23],1,1);lmiterm([-23 1 1 0],-R2);
% % pos (1, 2)
% lmiterm([-23 1 2 P23], A23, 1);
% % pos(2, 2)
% lmiterm([-23 2 2 P23], 1, 1);
% 
%  % pos in (1, 1)
%  lmiterm([-24 1 1 P24],1,1);lmiterm([-24 1 1 0],-R2);
% % % pos (1, 2)
%  lmiterm([-24 1 2 P24], A24, 1);
% % % pos(2, 2)
%  lmiterm([-24 2 2 P24], 1, 1);



