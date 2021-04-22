PA2=P*A2;
PA4=P2*A4;
PA6=P3*A6;
PA8=P4*A8;
PA10=P5*A10;
PA12=P6*A12;
PA14=P7*A14;
PA16=P8*A16;
PA18=P9*A18;
PA20=P10*A20;
PA22=P11*A22;
PA24=P12*A24;
PA26=P13*A26;
PA28=P14*A28;
PA30=P15*A30;
PA32=P16*A32;
PA34=P17*A34;
PA36=P18*A36;
PA38=P19*A38;
% PA40=P20*A40;
% PA42=P21*A42;
% PA44=P22*A44;
% PA46=P23*A46;
% PA48=P24*A48;

PN=-P;
E2=-eye(2);

% B22=-0.5858*B;%前面的系数，是L矩阵的特征值 B22=-B2;B44=-B4;B66=-B6,重复的，用处不大
% B44=-2*B;%
% B66=-3.4142*B;%

% % setlmis([]);  % DXW move to linesearch main
% p= lmivar(1, [2 1]);p2= lmivar(1, [2 1]);p3= lmivar(1, [2 1]);p4= lmivar(1, [2 1]);p5= lmivar(1, [2 1]);p6= lmivar(1, [2 1]);
% p7= lmivar(1, [2 1]);p8= lmivar(1, [2 1]);p9= lmivar(1, [2 1]);p10= lmivar(1, [2 1]);p11= lmivar(1, [2 1]);
% p12= lmivar(1, [2 1]);
% p13= lmivar(1, [2 1]);p14= lmivar(1, [2 1]);p15= lmivar(1, [2 1]);p16= lmivar(1, [2 1]);p17= lmivar(1, [2 1]);p18= lmivar(1, [2 1]);
% p19= lmivar(1, [2 1]);
% % p20= lmivar(1, [2 1]);p21= lmivar(1, [2 1]);p22= lmivar(1, [2 1]);p23= lmivar(1, [2 1]);p24= lmivar(1, [2 1]);
%  %   k= lmivar(2, [2 2]);  % DXW: moved to linesearch20shu main 
%描述矩阵p,P正定？
lmiterm([-4 1 1 p],1,1);
lmiterm([-4 1 1 0],P);
lmiterm([2 1 1 p],1,1);
lmiterm([2 1 1 0],E2);
lmiterm([-3 1 1 p],1,1);
lmiterm([-3 1 1 0],eye(2));

lmiterm([-6 1 1 p2],1,1);
lmiterm([-6 1 1 0],P2);
lmiterm([-5 1 1 p2],1,1);
lmiterm([-5 1 1 0],eye(2));

lmiterm([-10 1 1 p3],1,1);
lmiterm([-10 1 1 0],P3);
lmiterm([-9 1 1 p3],1,1);
lmiterm([-9 1 1 0],eye(2));

lmiterm([-12 1 1 p4],1,1);
lmiterm([-12 1 1 0],P4);
lmiterm([-11 1 1 p4],1,1);
lmiterm([-11 1 1 0],eye(2));

lmiterm([-15 1 1 p5],1,1);
lmiterm([-15 1 1 0],P5);
lmiterm([-14 1 1 p5],1,1);
lmiterm([-14 1 1 0],eye(2));

lmiterm([-18 1 1 p6],1,1);
lmiterm([-18 1 1 0],P6);
lmiterm([-17 1 1 p6],1,1);
lmiterm([-17 1 1 0],eye(2));

lmiterm([-21 1 1 p7],1,1);
lmiterm([-21 1 1 0],P7);
lmiterm([-20 1 1 p7],1,1);
lmiterm([-20 1 1 0],eye(2));

lmiterm([-24 1 1 p8],1,1);
lmiterm([-24 1 1 0],P8);
lmiterm([-23 1 1 p8],1,1);
lmiterm([-23 1 1 0],eye(2));

lmiterm([-27 1 1 p9],1,1);
lmiterm([-27 1 1 0],P9);
lmiterm([-26 1 1 p9],1,1);
lmiterm([-26 1 1 0],eye(2));

lmiterm([-30 1 1 p10],1,1);
lmiterm([-30 1 1 0],P10);
lmiterm([-29 1 1 p10],1,1);
lmiterm([-29 1 1 0],eye(2));

lmiterm([-33 1 1 p11],1,1);
lmiterm([-33 1 1 0],P11);
lmiterm([-32 1 1 p11],1,1);
lmiterm([-32 1 1 0],eye(2));

lmiterm([-36 1 1 p12],1,1);
lmiterm([-36 1 1 0],P12);
lmiterm([-36 1 1 p12],1,1);
lmiterm([-36 1 1 0],eye(2));

lmiterm([-39 1 1 p13],1,1);
lmiterm([-39 1 1 0],P13);
lmiterm([-38 1 1 p13],1,1);
lmiterm([-38 1 1 0],eye(2));

lmiterm([-42 1 1 p14],1,1);
lmiterm([-42 1 1 0],P14);
lmiterm([-41 1 1 p14],1,1);
lmiterm([-41 1 1 0],eye(2));

lmiterm([-45 1 1 p15],1,1);
lmiterm([-45 1 1 0],P15);
lmiterm([-44 1 1 p15],1,1);
lmiterm([-44 1 1 0],eye(2));

lmiterm([-48 1 1 p16],1,1);
lmiterm([-48 1 1 0],P16);
lmiterm([-47 1 1 p16],1,1);
lmiterm([-47 1 1 0],eye(2));

lmiterm([-51 1 1 p17],1,1);
lmiterm([-51 1 1 0],P17);
lmiterm([-50 1 1 p17],1,1);
lmiterm([-50 1 1 0],eye(2));

lmiterm([-54 1 1 p18],1,1);
lmiterm([-54 1 1 0],P18);
lmiterm([-53 1 1 p18],1,1);
lmiterm([-53 1 1 0],eye(2));

lmiterm([-57 1 1 p19],1,1);
lmiterm([-57 1 1 0],P19);
lmiterm([-56 1 1 p19],1,1);
lmiterm([-56 1 1 0],eye(2));
% 
% lmiterm([-60 1 1 p20],1,1);
% lmiterm([-60 1 1 0],P20);
% lmiterm([-59 1 1 p20],1,1);
% lmiterm([-59 1 1 0],eye(2));
% 
% lmiterm([-63 1 1 p21],1,1);
% lmiterm([-63 1 1 0],P21);
% lmiterm([-62 1 1 p21],1,1);
% lmiterm([-62 1 1 0],eye(2));
% 
% lmiterm([-66 1 1 p22],1,1);
% lmiterm([-66 1 1 0],P22);
% lmiterm([-65 1 1 p22],1,1);
% lmiterm([-65 1 1 0],eye(2));
% 
% lmiterm([-69 1 1 p23],1,1);
% lmiterm([-69 1 1 0],P23);
% lmiterm([-68 1 1 p23],1,1);
% lmiterm([-68 1 1 0],eye(2));
% 
% lmiterm([-72 1 1 p24],1,1);
% lmiterm([-72 1 1 0],P24);
% lmiterm([-71 1 1 p24],1,1);
% lmiterm([-71 1 1 0],eye(2));






lmiterm([-1 1 1 p],1,1);%p-A1*p*A2+B2*k*PA2+PA2'*k'*B2'
lmiterm([-1 1 1 p],A1,-A2);
lmiterm([-1 1 1 k],B2,PA2,'s');
lmiterm([-1 1 2 p],A1,1);%A1*p-B*lamda2*k
lmiterm([-1 1 2 k],-B,lamda2);
lmiterm([-1 1 3 k],B2,P);%B2*k*P=lamda2*B*k*P
lmiterm([-1 2 2 0],eye(2));%I
lmiterm([-1 3 3 0],P);%P

lmiterm([-2 1 1 p2],1,1);
lmiterm([-2 1 1 p2],A3,-A4);
lmiterm([-2 1 1 k],B4,PA4,'s');
lmiterm([-2 1 2 p2],A3,1);
lmiterm([-2 1 2 k],-B,lamda3);
lmiterm([-2 1 3 k],B4,P2);
lmiterm([-2 2 2 0],eye(2));
lmiterm([-2 3 3 0],P2);

lmiterm([-3 1 1 p3],1,1);
lmiterm([-3 1 1 p3],A5,-A6);
lmiterm([-3 1 1 k],B6,PA6,'s');
lmiterm([-3 1 2 p3],A5,1);
lmiterm([-3 1 2 k],-B,lamda4);
lmiterm([-3 1 3 k],B6,P3);
lmiterm([-3 2 2 0],eye(2));
lmiterm([-3 3 3 0],P3);

lmiterm([-4 1 1 p4],1,1);
lmiterm([-4 1 1 p4],A7,-A8);
lmiterm([-4 1 1 k],B8,PA8,'s');
lmiterm([-4 1 2 p4],A7,1);
lmiterm([-4 1 2 k],-B,lamda5);
lmiterm([-4 1 3 k],B8,P4);
lmiterm([-4 2 2 0],eye(2));
lmiterm([-4 3 3 0],P4);

lmiterm([-5 1 1 p5],1,1);
lmiterm([-5 1 1 p5],A9,-A10);
lmiterm([-5 1 1 k],B10,PA10,'s');
lmiterm([-5 1 2 p5],A9,1);
lmiterm([-5 1 2 k],-B,lamda6);
lmiterm([-5 1 3 k],B10,P5);
lmiterm([-5 2 2 0],eye(2));
lmiterm([-5 3 3 0],P5);

lmiterm([-6 1 1 p6],1,1);
lmiterm([-6 1 1 p6],A11,-A12);
lmiterm([-6 1 1 k],B12,PA12,'s');
lmiterm([-6 1 2 p6],A11,1);
lmiterm([-6 1 2 k],-B,lamda7);
lmiterm([-6 1 3 k],B12,P6);
lmiterm([-6 2 2 0],eye(2));
lmiterm([-6 3 3 0],P6);

lmiterm([-7 1 1 p7],1,1);
lmiterm([-7 1 1 p7],A13,-A14);
lmiterm([-7 1 1 k],B14,PA14,'s');
lmiterm([-7 1 2 p7],A13,1);
lmiterm([-7 1 2 k],-B,lamda8);
lmiterm([-7 1 3 k],B14,P7);
lmiterm([-7 2 2 0],eye(2));
lmiterm([-7 3 3 0],P7);

lmiterm([-8 1 1 p8],1,1);
lmiterm([-8 1 1 p8],A15,-A16);
lmiterm([-8 1 1 k],B16,PA16,'s');
lmiterm([-8 1 2 p8],A15,1);
lmiterm([-8 1 2 k],-B,lamda9);
lmiterm([-8 1 3 k],B16,P8);
lmiterm([-8 2 2 0],eye(2));
lmiterm([-8 3 3 0],P8);

lmiterm([-9 1 1 p9],1,1);
lmiterm([-9 1 1 p9],A17,-A18);
lmiterm([-9 1 1 k],B18,PA18,'s');
lmiterm([-9 1 2 p9],A17,1);
lmiterm([-9 1 2 k],-B,lamda10);
lmiterm([-9 1 3 k],B18,P9);
lmiterm([-9 2 2 0],eye(2));
lmiterm([-9 3 3 0],P9);

lmiterm([-10 1 1 p10],1,1);
lmiterm([-10 1 1 p10],A19,-A20);
lmiterm([-10 1 1 k],B20,PA20,'s');
lmiterm([-10 1 2 p10],A19,1);
lmiterm([-10 1 2 k],-B,lamda11);
lmiterm([-10 1 3 k],B20,P10);
lmiterm([-10 2 2 0],eye(2));
lmiterm([-10 3 3 0],P10);

lmiterm([-11 1 1 p11],1,1);
lmiterm([-11 1 1 p11],A21,-A22);
lmiterm([-11 1 1 k],B22,PA22,'s');
lmiterm([-11 1 2 p11],A21,1);
lmiterm([-11 1 2 k],-B,lamda12);
lmiterm([-11 1 3 k],B22,P11);
lmiterm([-11 2 2 0],eye(2));
lmiterm([-11 3 3 0],P11);

lmiterm([-12 1 1 p12],1,1);
lmiterm([-12 1 1 p12],A23,-A24);
lmiterm([-12 1 1 k],B24,PA24,'s');
lmiterm([-12 1 2 p12],A23,1);
lmiterm([-12 1 2 k],-B,lamda13);
lmiterm([-12 1 3 k],B24,P12);
lmiterm([-12 2 2 0],eye(2));
lmiterm([-12 3 3 0],P12);

lmiterm([-13 1 1 p13],1,1);
lmiterm([-13 1 1 p13],A25,-A26);
lmiterm([-13 1 1 k],B26,PA26,'s');
lmiterm([-13 1 2 p13],A25,1);
lmiterm([-13 1 2 k],-B,lamda14);
lmiterm([-13 1 3 k],B26,P13);
lmiterm([-13 2 2 0],eye(2));
lmiterm([-13 3 3 0],P13);

lmiterm([-14 1 1 p14],1,1);
lmiterm([-14 1 1 p14],A27,-A28);
lmiterm([-14 1 1 k],B28,PA28,'s');
lmiterm([-14 1 2 p14],A27,1);
lmiterm([-14 1 2 k],-B,lamda15);
lmiterm([-14 1 3 k],B28,P14);
lmiterm([-14 2 2 0],eye(2));
lmiterm([-14 3 3 0],P14);

lmiterm([-15 1 1 p15],1,1);
lmiterm([-15 1 1 p15],A29,-A30);
lmiterm([-15 1 1 k],B30,PA30,'s');
lmiterm([-15 1 2 p15],A29,1);
lmiterm([-15 1 2 k],-B,lamda16);
lmiterm([-15 1 3 k],B30,P15);
lmiterm([-15 2 2 0],eye(2));
lmiterm([-15 3 3 0],P15);

lmiterm([-16 1 1 p16],1,1);
lmiterm([-16 1 1 p16],A31,-A32);
lmiterm([-16 1 1 k],B32,PA32,'s');
lmiterm([-16 1 2 p16],A31,1);
lmiterm([-16 1 2 k],-B,lamda17);
lmiterm([-16 1 3 k],B32,P16);
lmiterm([-16 2 2 0],eye(2));
lmiterm([-16 3 3 0],P16);

lmiterm([-17 1 1 p17],1,1);
lmiterm([-17 1 1 p17],A33,-A34);
lmiterm([-17 1 1 k],B34,PA34,'s');
lmiterm([-17 1 2 p17],A33,1);
lmiterm([-17 1 2 k],-B,lamda18);
lmiterm([-17 1 3 k],B34,P17);
lmiterm([-17 2 2 0],eye(2));
lmiterm([-17 3 3 0],P17);

lmiterm([-18 1 1 p18],1,1);
lmiterm([-18 1 1 p18],A35,-A36);
lmiterm([-18 1 1 k],B36,PA36,'s');
lmiterm([-18 1 2 p18],A35,1);
lmiterm([-18 1 2 k],-B,lamda19);
lmiterm([-18 1 3 k],B36,P18);
lmiterm([-18 2 2 0],eye(2));
lmiterm([-18 3 3 0],P18);

lmiterm([-19 1 1 p19],1,1);
lmiterm([-19 1 1 p19],A37,-A38);
lmiterm([-19 1 1 k],B38,PA38,'s');
lmiterm([-19 1 2 p19],A37,1);
lmiterm([-19 1 2 k],-B,lamda20);
lmiterm([-19 1 3 k],B38,P19);
lmiterm([-19 2 2 0],eye(2));
lmiterm([-19 3 3 0],P19);
% 
% lmiterm([-20 1 1 p20],1,1);
% lmiterm([-20 1 1 p20],A39,A40);
% lmiterm([-20 1 1 k],B40,PA40,'s');
% lmiterm([-20 1 2 p20],A39,1);
% lmiterm([-20 1 2 k],-B,2);
% lmiterm([-20 1 3 k],B39,P20);
% lmiterm([-20 2 2 0],eye(2));
% lmiterm([-20 3 3 0],P20);
% 
% lmiterm([-21 1 1 p21],1,1);
% lmiterm([-21 1 1 p21],A41,A42);
% lmiterm([-21 1 1 k],B42,PA42,'s');
% lmiterm([-21 1 2 p21],A41,1);
% lmiterm([-21 1 2 k],-B,2);
% lmiterm([-21 1 3 k],B41,P21);
% lmiterm([-21 2 2 0],eye(2));
% lmiterm([-21 3 3 0],P21);
% 
% lmiterm([-22 1 1 p22],1,1);
% lmiterm([-22 1 1 p22],A43,A44);
% lmiterm([-22 1 1 k],B44,PA44,'s');
% lmiterm([-22 1 2 p22],A43,1);
% lmiterm([-22 1 2 k],-B,2);
% lmiterm([-22 1 3 k],B43,P22);
% lmiterm([-22 2 2 0],eye(2));
% lmiterm([-22 3 3 0],P22);
% 
% lmiterm([-23 1 1 p23],1,1);
% lmiterm([-23 1 1 p23],A45,A46);
% lmiterm([-23 1 1 k],B46,PA46,'s');
% lmiterm([-23 1 2 p23],A45,1);
% lmiterm([-23 1 2 k],-B,2);
% lmiterm([-23 1 3 k],B45,P23);
% lmiterm([-23 2 2 0],eye(2));
% lmiterm([-23 3 3 0],P23);
% 
% lmiterm([-24 1 1 p24],1,1);
% lmiterm([-24 1 1 p24],A47,A48);
% lmiterm([-24 1 1 k],B48,PA48,'s');
% lmiterm([-24 1 2 p24],A47,1);
% lmiterm([-24 1 2 k],-B,2);
% lmiterm([-24 1 3 k],B47,P24);
% lmiterm([-24 2 2 0],eye(2));
% lmiterm([-24 3 3 0],P24);