        Akasi=A-lamda2*B*(K+kasi1*deltak);
        Akasi2=A-lamda3*B*(K+kasi1*deltak);
        Akasi3=A-lamda4*B*(K+kasi1*deltak);
        Akasi4=A-lamda5*B*(K+kasi1*deltak);
        Akasi5=A-lamda6*B*(K+kasi1*deltak);
        Akasi6=A-lamda7*B*(K+kasi1*deltak);
        Akasi7=A-lamda8*B*(K+kasi1*deltak);
        Akasi8=A-lamda9*B*(K+kasi1*deltak);
        Akasi9=A-lamda10*B*(K+kasi1*deltak);
        Akasi10=A-lamda11*B*(K+kasi1*deltak);
        Akasi11=A-lamda12*B*(K+kasi1*deltak);
        Akasi12=A-lamda13*B*(K+kasi1*deltak);
        Akasi13=A-lamda14*B*(K+kasi1*deltak);
        Akasi14=A-lamda15*B*(K+kasi1*deltak);
        Akasi15=A-lamda16*B*(K+kasi1*deltak);
        Akasi16=A-lamda17*B*(K+kasi1*deltak);
        Akasi17=A-lamda18*B*(K+kasi1*deltak);
        Akasi18=A-lamda19*B*(K+kasi1*deltak);
        Akasi19=A-lamda20*B*(K+kasi1*deltak);
        
        
        
        %     Akasi1=-Akasi;
        setlmis([]);
        P= lmivar(1, [2 1]);P2= lmivar(1, [2 1]);P3= lmivar(1, [2 1]);P4 = lmivar(1, [2 1]);P5 = lmivar(1, [2 1]);P6 = lmivar(1, [2 1]);
        P7 = lmivar(1, [2 1]);P8 = lmivar(1, [2 1]);P9 = lmivar(1, [2 1]);P10 = lmivar(1, [2 1]);P11 = lmivar(1, [2 1]);P12 = lmivar(1, [2 1]);
        P13 = lmivar(1, [2 1]);P14 = lmivar(1, [2 1]);P15 = lmivar(1, [2 1]);P16 = lmivar(1, [2 1]);P17 = lmivar(1, [2 1]);P18 = lmivar(1, [2 1]);
        P19 = lmivar(1, [2 1]);
        %      P20 = lmivar(1, [2 1]);P21 = lmivar(1, [2 1]);P22 = lmivar(1, [2 1]);P23 = lmivar(1, [2 1]);P24 = lmivar(1, [2 1]);
        
        %     lmiterm([-2 1 1 P],1,1);
        %     lmiterm([-5 1 1 P2],1,1);
        %     lmiterm([-6 1 1 P3],1,1);
        lmiterm([-2 1 1 P], 1, 1);lmiterm([-5 1 1 P2], 1, 1);lmiterm([-6 1 1 P3], 1, 1);lmiterm([-2 1 1 P4], 1, 1);lmiterm([-5 1 1 P5], 1, 1);lmiterm([-6 1 1 P6], 1, 1);
        lmiterm([-2 1 1 P7], 1, 1);lmiterm([-5 1 1 P8], 1, 1);lmiterm([-6 1 1 P9], 1, 1);lmiterm([-2 1 1 P10], 1, 1);lmiterm([-5 1 1 P11], 1, 1);lmiterm([-6 1 1 P12], 1, 1);
        lmiterm([-2 1 1 P13], 1, 1);lmiterm([-5 1 1 P14], 1, 1);lmiterm([-6 1 1 P15], 1, 1);lmiterm([-2 1 1 P16], 1, 1);lmiterm([-5 1 1 P17], 1, 1);lmiterm([-6 1 1 P18], 1, 1);
        lmiterm([-2 1 1 P19], 1, 1);
        %     lmiterm([-5 1 1 P20], 1, 1);lmiterm([-6 1 1 P21], 1, 1);lmiterm([-2 1 1 P22], 1, 1);lmiterm([-5 1 1 P23], 1, 1);lmiterm([-6 1 1 P24], 1, 1);
        
        %     lmiterm([-1 1 1 P],1,1);
        %     lmiterm([-1 1 1 0],-R2);
        %     lmiterm([-1 1 2 P], Akasi, 1);
        %     lmiterm([-1 2 2 P], 1, 1);
        %
        %     lmiterm([-3 1 1 P2],1,1);
        %     lmiterm([-3 1 1 0],-R2);
        %     lmiterm([-3 1 2 P2], Akasi2, 1);
        %     lmiterm([-3 2 2 P2], 1, 1);
        %
        %     lmiterm([-4 1 1 P2],1,1);
        %     lmiterm([-4 1 1 0],-R2);
        %     lmiterm([-4 1 2 P2], Akasi3, 1);
        %     lmiterm([-4 2 2 P2], 1, 1);
        
        
        lmiterm([-1 1 1 P],1,1);lmiterm([-1 1 1 0],-R2);
        lmiterm([-1 1 2 P], Akasi, 1);
        lmiterm([-1 2 2 P], 1, 1);
        
        lmiterm([-2 1 1 P2],1,1);lmiterm([-2 1 1 0],-R2);
        lmiterm([-2 1 2 P2], Akasi2, 1);
        lmiterm([-2 2 2 P2], 1, 1);
        
        lmiterm([-3 1 1 P3],1,1);lmiterm([-3 1 1 0],-R2);
        lmiterm([-3 1 2 P3], Akasi3, 1);
        lmiterm([-3 2 2 P3], 1, 1);
        
        lmiterm([-4 1 1 P4],1,1);lmiterm([-4 1 1 0],-R2);
        lmiterm([-4 1 2 P4], Akasi4, 1);
        lmiterm([-4 2 2 P4], 1, 1);
        
        lmiterm([-5 1 1 P5],1,1);lmiterm([-5 1 1 0],-R2);
        lmiterm([-5 1 2 P5],Akasi5, 1);
        lmiterm([-5 2 2 P5], 1, 1);
        
        lmiterm([-6 1 1 P6],1,1);lmiterm([-6 1 1 0],-R2);
        lmiterm([-6 1 2 P6], Akasi6, 1);
        lmiterm([-6 2 2 P6], 1, 1);
        
        
        lmiterm([-7 1 1 P7],1,1);lmiterm([-7 1 1 0],-R2);
        lmiterm([-7 1 2 P7], Akasi7, 1);
        lmiterm([-7 2 2 P7], 1, 1);
        
        lmiterm([-8 1 1 P8],1,1);lmiterm([-8 1 1 0],-R2);
        lmiterm([-8 1 2 P8], Akasi8, 1);
        lmiterm([-8 2 2 P8], 1, 1);
        
        lmiterm([-9 1 1 P9],1,1);lmiterm([-3 1 1 0],-R2);
        lmiterm([-9 1 2 P9], Akasi9, 1);
        lmiterm([-9 2 2 P9], 1, 1);
        
        lmiterm([-10 1 1 P10],1,1);lmiterm([-10 1 1 0],-R2);
        lmiterm([-10 1 2 P10], Akasi10, 1);
        lmiterm([-10 2 2 P10], 1, 1);
        
        lmiterm([-11 1 1 P11],1,1);lmiterm([-11 1 1 0],-R2);
        lmiterm([-11 1 2 P11], Akasi11, 1);
        lmiterm([-11 2 2 P11], 1, 1);
        
        lmiterm([-12 1 1 P12],1,1);lmiterm([-12 1 1 0],-R2);
        lmiterm([-12 1 2 P12], Akasi12, 1);
        lmiterm([-12 2 2 P12], 1, 1);
        
        lmiterm([-13 1 1 P13],1,1);lmiterm([-13 1 1 0],-R2);
        lmiterm([-13 1 2 P13], Akasi13, 1);
        lmiterm([-13 2 2 P13], 1, 1);
        
        lmiterm([-14 1 1 P14],1,1);lmiterm([-14 1 1 0],-R2);
        lmiterm([-14 1 2 P14], Akasi14, 1);
        lmiterm([-14 2 2 P14], 1, 1);
        
        lmiterm([-15 1 1 P15],1,1);lmiterm([-15 1 1 0],-R2);
        lmiterm([-15 1 2 P15], Akasi15, 1);
        lmiterm([-15 2 2 P15], 1, 1);
        
        lmiterm([-16 1 1 P16],1,1);lmiterm([-16 1 1 0],-R2);
        lmiterm([-16 1 2 P16], Akasi16, 1);
        lmiterm([-16 2 2 P16], 1, 1);
        
        lmiterm([-17 1 1 P17],1,1);lmiterm([-17 1 1 0],-R2);
        lmiterm([-17 1 2 P17], Akasi17, 1);
        lmiterm([-17 2 2 P17], 1, 1);
        
        lmiterm([-18 1 1 P18],1,1);lmiterm([-18 1 1 0],-R2);
        lmiterm([-18 1 2 P18], Akasi18, 1);
        lmiterm([-18 2 2 P18], 1, 1);
        
        lmiterm([-19 1 1 P19],1,1);lmiterm([-19 1 1 0],-R2);
        lmiterm([-19 1 2 P19], Akasi19, 1);
        lmiterm([-19 2 2 P19], 1, 1);
        %
        % lmiterm([-20 1 1 P20],1,1);lmiterm([-20 1 1 0],-R2);
        % lmiterm([-20 1 2 P20], Akasi20, 1);
        % lmiterm([-20 2 2 P20], 1, 1);
        %
        % lmiterm([-21 1 1 P21],1,1);lmiterm([-21 1 1 0],-R2);
        % lmiterm([-21 1 2 P21], Akasi21, 1);
        % lmiterm([-21 2 2 P21], 1, 1);
        %
        % lmiterm([-22 1 1 P22],1,1);lmiterm([-22 1 1 0],-R2);
        % lmiterm([-22 1 2 P22], Akasi22, 1);
        % lmiterm([-22 2 2 P22], 1, 1);
        %
        % lmiterm([-23 1 1 P23],1,1);lmiterm([-23 1 1 0],-R2);
        % lmiterm([-23 1 2 P23], Akasi23, 1);
        % lmiterm([-23 2 2 P23], 1, 1);
        %
        % % lmiterm([-24 1 1 P24],1,1);lmiterm([-24 1 1 0],-R2);
        % % lmiterm([-24 1 2 P24], Akasi24, 1);
        % % lmiterm([-24 2 2 P24], 1, 1);
        
        