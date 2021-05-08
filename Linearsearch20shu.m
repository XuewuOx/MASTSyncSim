% clear all
% clc
if ~exist('DEBUGPLOT','var')
    DEBUGPLOT = false; % switch to plot medium results for debug
end

disp(" ****  Opimize K for Mimimum Variance  **** ");
if ~exist('A','var') ||~exist('B','var')
    A=[1 1;0 1];
    B=[1 1;0 1];
end
% set the Laplacian matrix L. On 0416,default L=Lshu20
if ~exist('L','var') || isempty(L)
    fprintf("L is empty. Create a new network L by genNetbyL()\n");
    [netG L]=genNetbyL(); 
end
%偶数是特征值*B，奇数是（-前一项）特征值：0；3.4142；2；0.5858；1（20个）
lamda1=0;lamda2=0.0698;lamda3=0.1390;lamda4=0.2679;lamda5=0.2679;lamda6=0.5967;lamda7=0.9344;lamda8=1.0000;lamda9=1.0000;lamda10=1.0000;lamda11=1.0000;lamda12=1.7853;lamda13=2.2057;
lamda14=2.6042;lamda15=3.3312;lamda16=3.7321;lamda17=3.7321;lamda18=4.3475;lamda19=4.8322;lamda20=5.1541;

B2=lamda2*B;B4=lamda3*B;B6=lamda4*B;B8=lamda5*B;B10=lamda6*B ;B12=lamda7*B;B14=lamda8*B;B16=lamda9*B;B18=lamda10*B;B20=lamda11*B;B22=lamda12*B;B24=lamda13*B;B26=lamda14*B;B28=lamda15*B;B30=lamda16*B;B32=lamda17*B;B34=lamda18*B;B36=lamda19*B;B38=lamda20*B;
B3=-B2;B5=-B4;B7=-B6;B9=-B8;B11=-B10;B13=-B12;B15=-B14;B17=-B16;B19=-B18;B21=-B20;B23=-B22;B25=-B24;B27=-B26;B29=-B28;B31=-B30;B33=-B32;B35=-B34;B37=-B36;B39=-B38;
% B2=1*B;B4=1*B;B6=1*B;B8=1*B;B10=1*B;B12=1*B;B14=1*B;B16=1*B;B18=1*B;B20=1*B;B22=1*B;B24=1*B;B26=1*B;B28=1*B;B30=1*B;B32=1*B;B34=1*B;B36=1*B;B38=1*B;B40=1*B;
% B3=-B2;B5=-B4;B7=-B6;B9=-B8;B11=-B10;B13=-B12;B15=-B14;B17=-B16;B19=-B18;B21=-B20;B23=-B22;B25=-B24;B27=-B26;B29=-B28;B31=-B30;B33=-B32;B35=-B34;B37=-B36;B39=-B38;B41=-B40;
% B42=3.4142*B;
% B43=-B2;
% B44=2*B;
% B45=-B4;
% B46=0.5858*B;
% B47=-B6;

% check if noise covariance matrix R Q have been defined
if ~exist('R','var') || ~exist('Q','var') 
  % Q, R not defined, set with default values
  R=10^-16*[1 0;0 1];%0416过程噪声
  Q=10^-14*[1 1;1 2];%测量噪声
end

fprintf("Process Noise Cov R="); disp(R);
fprintf("Observation Noise Cov Q ="); disp(Q);

vecOne=ones(4,1);
l=eye(4)-1/4*(vecOne*vecOne');
% U=[-0.6533 -0.500 -0.2706;-0.2706 0.500 0.6533;0.2706 0.500 -0.6533;0.6533 -0.5000 0.2706];
R1=A*Q*A';
R2=R+R1+Q;

% select the initial value of K
if exist('K0','var')
    % if K0 is already given
    K=K0;
else
    % ninitial value of K is not specified, 
    K01=[0.1 0;0 0.1];%原来的初始值
    K02=[0.3517 -0.1997;0.0159 0.2565];%以传统中得到的优化后的增益作为初始值,J变大了10-14
    K03=[0.3338 -0.0840;0.0133 0.3034];%以传统中得到的优化后的增益作为初始值 10-14  10-16

    K=K01;
end
 disp([' * Initial value of K is ']); disp(K);
 
 disp(['   Check the stability of gain K against A,B,L']);
[stableFlag,eigAnetc,numEigOne] = chkEigAc(A,B,K,L);
if (stableFlag==false)
    warning("Non-stable control gain matrix");
    % return; % debug codes
    nID=input('Would you like to contintue with the Non-stable gain matrix? \n','s');
    if ~(nID=='y' ||  nID=='Y')
        return;
    end
end

taoe=10^-10;
size_kchangex=1;%定义的一个中间变量
taok=10^-10;
szopt=35; % max interation of main optimization steps
szstepsize=20; % max interation of stepsize optimization


% J=zeros(1,sz1);
%%
% % %A*的奇数行是A-lamdaBK,A*的偶数行是A-1的转置
% % A1=A-B2*K;A2=A1';
% % A3=A-B4*K;A4=A3';
% % A5=A-B6*K;A6=A5';
% % A7=A-B8*K;A8=A7';
% % A9=A-B10*K;A10=A9';
% % A11=A-B12*K;A12=A11';
% % A13=A-B14*K;A14=A13';
% % A15=A-B16*K;A16=A15';
% % A17=A-B18*K;A18=A17';
% % A19=A-B20*K;A20=A19';
% % A21=A-B22*K;A22=A21';
% % A23=A-B24*K;A24=A23';
% % A25=A-B26*K;A26=A25';
% % A27=A-B28*K;A28=A27';
% % A29=A-B30*K;A30=A29';
% % A31=A-B32*K;A32=A31';
% % A33=A-B34*K;A34=A33';
% % A35=A-B36*K;A36=A35';
% % A37=A-B38*K;A38=A37';
% % % A39=A-B40*K;A40=A39';
% % % A41=A-B42*K;A42=A41';
% % % A43=A-B*K;A44=A43';
% % % A45=A-B*K;A46=A45';
% % % A47=A-B*K;A48=A47';
% % 

%%
if DEBUGPLOT
    hfig_Popt=figure('name','Optimizing tr(P1+P2+...+PN-1)');
    xlabel("main iteration i"); ylabel("Ycopt"); hold on;
end

J1x=[];
allkasi=[]; % to save the value of allkasi
allk2=[];
allK=[];
% sz=10;
for i=1:szopt
     
    fprintf(" *  i=%d \n",i);
    % solve eq.(40) find optimal tr(P1+P2+...+PN-1)
    disp("   Step 1: Solve eq.(40) to find optimal tr(P1+P2+...+PN-1)");
    setLMIProb_Popt;
    
    lmis = getlmis;
    c=mat2dec(lmis,eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2));
    options=[1e-8 0 0 0 1]; % options(5)=0 turn on execution trace; =1 turn off execution trace
    [copt, xopt]=mincx(lmis,c,options);
    fprintf("      best objective value: J1(%d)= %d \n", i, copt);
    J1=copt;  J1x(i)=copt;
    if DEBUGPLOT
        figure(hfig_Popt); plot(J1x,'-ob','LineWidth',2);
    end
    P=dec2mat(lmis,xopt,P);P2=dec2mat(lmis,xopt,P2);P3=dec2mat(lmis,xopt,P3);P4=dec2mat(lmis,xopt,P4);P5=dec2mat(lmis,xopt,P5);P6=dec2mat(lmis,xopt,P6);
    P7=dec2mat(lmis,xopt,P7);P8=dec2mat(lmis,xopt,P8);P9=dec2mat(lmis,xopt,P9);P10=dec2mat(lmis,xopt,P10);P11=dec2mat(lmis,xopt,P11);
    P12=dec2mat(lmis,xopt,P12);
    P13=dec2mat(lmis,xopt,P13);P14=dec2mat(lmis,xopt,P14);P15=dec2mat(lmis,xopt,P15);P16=dec2mat(lmis,xopt,P16);P17=dec2mat(lmis,xopt,P17);P18=dec2mat(lmis,xopt,P18);
    P19=dec2mat(lmis,xopt,P19);%P20=dec2mat(lmis,xopt,P20);P21=dec2mat(lmis,xopt,P21);P22=dec2mat(lmis,xopt,P22);P23=dec2mat(lmis,xopt,P23);P24=dec2mat(lmis,xopt,P24);
    % J(1)=P(1,1)+P(2,2); % for 2 nodes 
   
    %%
    % solve eq.(44) find Delta K to minimize  tr(Delta P1 + Delta P2+...+ Delta PN-1)
    disp("   Step 2: Solve eq.(44) to find optimal tr(Delta P1+ Delta P2+...+ Delta PN-1)");
    setlmis([]);
    p= lmivar(1, [2 1]);p2= lmivar(1, [2 1]);p3= lmivar(1, [2 1]);p4= lmivar(1, [2 1]);p5= lmivar(1, [2 1]);p6= lmivar(1, [2 1]);
    p7= lmivar(1, [2 1]);p8= lmivar(1, [2 1]);p9= lmivar(1, [2 1]);p10= lmivar(1, [2 1]);p11= lmivar(1, [2 1]);
    p12= lmivar(1, [2 1]);
    p13= lmivar(1, [2 1]);p14= lmivar(1, [2 1]);p15= lmivar(1, [2 1]);p16= lmivar(1, [2 1]);p17= lmivar(1, [2 1]);p18= lmivar(1, [2 1]);
    p19= lmivar(1, [2 1]);
    % p20= lmivar(1, [2 1]);p21= lmivar(1, [2 1]);p22= lmivar(1, [2 1]);p23= lmivar(1, [2 1]);p24= lmivar(1, [2 1]);
    k= lmivar(2, [2 2]); % k is used in InitLMIProb_deltaK
    setLMIProb_deltaK; % rename to InitLMIProb_deltaK
    
    lmis1 = getlmis;
    c1=mat2dec(lmis1,eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),zeros(2,2));
    options1=[1e-8 0 0 0 1];
    [copt1,xopt1]=mincx(lmis1,c1,options1);
    if isempty(copt1)
         warning("      no feasible solution of deltaK found. Continue to next episode without change K\n");
         continue;
    end
    deltak=dec2mat(lmis1,xopt1,k);
    fprintf("      best objective value deltaK: %d at deltaK=\n",copt1);
    if DEBUGPLOT  disp(deltak); end;
    
    
    % K=K+ Kasi* Delta K
    
    %% Step 3: Linear Search for optimal stepsize
    disp("   Step 3: Linear Search for optimal stepsize");
    % J=zeros(1,sz1);
    %     Akasi=A-2*B*(K+kasi*k);
    %     Akasi1=-Akasi;
    %     setlmis([]);
    %     P= lmivar(1, [2 1]);
    %     lmiterm([-2 1 1 P],1,1);
    %     lmiterm([1 1 1 P],-1,1);
    %     lmiterm([1 1 1 0],R2);
    %     lmiterm([1 1 2 P], Akasi1, 1);
    %     lmiterm([1 2 2 P], -1, 1);
    %     lmis3 = getlmis;
    %     c3=mat2dec(lmis3,eye(2));
    %     options3=[1e-8 0 0 100 1];
    %     xinit3=zeros(3,1);
    %     [copt3 xopt3]=mincx(lmis3,c3,options3,xinit3);
    %     J(1)=xopt3(1)+xopt3(3);
    % if k2<taok
    %     break;
    % else
    %      kasi0=0.01;kasi=0;
    %      kasi=kasi0+kasi;
    % end
    kasi0=0.05;kasi=0; kasi1=0;
    kasix=[]; kasi0x=[];kasi1x=[];
    J2x=[]; % J1;
    Jminx=[]; % J1;
    Jmin=J1;
    szstepsize=13;
    % K=K*0.9;
    for j=1:szstepsize
       % plot(j,kasi,'ob'); plot(j,kasi0,'*r');
        kasix=[kasix kasi]; 
       kasi1x=[kasi1x kasi1];
               
        setLMIProb_lineStep
        
        lmis4 = getlmis;
        c4=mat2dec(lmis4,eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2),eye(2));
        options4=[1e-8 0 0 0 1];
        [copt4, xopt4]=mincx(lmis4,c4,options4);
        %     P=dec2mat(lmis4,xopt4,P);
        %     J(j)=P(1,1)+P(2,2);
        % TODO: check if a feasible solution is found. if not, set a failure flag
        % and set kasi to the last feasible value or the initial value 0
        if isempty(copt4)
             warning("     stepsizeOpt (j=%d) no feasible solution of line step found at Kasi1=%d\n", j,kasi1);
            fprintf("         Reduce the stepsize Kasi0 (%d) and try again \n", kasi0);
            J2=inf;  % set J2 to inf so that Ksai will not be updated
            J2x=[J2x nan]; % push nan (not inf) for better visualization 
        
%              kasi0x=[kasi0x kasi0];
%             Jminx=[Jminx Jmin];
%             kasi1=kasi+kasi0;
%             continue;
        else
            J2=copt4;
            J2x=[J2x J2];
        end
       
        if abs(kasi0)<taoe%步长优化
            break;
        else
            if J2<=Jmin
                kasi=kasi1;    % update kasi with new value
                % kasi=kasi1;%kasi1=kasi0+kasi;
                kasi0=2*kasi0;
                Jmin=J2;
            else
                % kasi=kasi+kasi0;    % update kasi with new value
                % kasi0=-0.25*kasi0;
                kasi0=-0.25*kasi0;
            end
        end
        kasi0x=[kasi0x kasi0];
        Jminx=[Jminx Jmin];
        kasi1=kasi+kasi0;
    end % for j=2:szstepsize , step size optimization
    
    J1=Jmin; % for optimize tr(P1+P2+...PN-1)
    if DEBUGPLOT==true
        figure('name',['Step size Kasi & J2 (i=' num2str(i) ')']);
        subplot(2,1,1);
        plot(kasix,'*-b');title(['Step size Kasix (blue) & kasi0x(red) (i=' num2str(i) ')']);
        hold on; plot(kasi0x,'o-r');
        % figure; plot(kasi0x,'*-r');title('Kasi0');
        % figure; plot(kasi1x,'o-k'); title('Kasi1');
        % figure('name', ['J2 in linear stepsize search (i=' num2str(i) ')']);
        subplot(2,1,2);
        plot(J2x,'^-m'); title(['J2 (pink) and Jmin (J1, black) in linear stepsize search (i=' num2str(i) ')']);
        hold on; plot(Jminx,'o-k'); 
    end
    
    %%  Update K
    %  save ('Step','kasi','-ascii','-append') %步长储存
    allkasi=[allkasi kasi];
    
    k1=abs(kasi*deltak);
    % k2=sqrt((k1(1,1)^2+k1(1,2)^2+k1(2,1)^2+k1(2,2)^2)/4);
    size_kchange=norm(k1);
    size_kchangex=[size_kchangex size_kchange];
    allK(i,:,:)=K;
    %  k2=sqrt(k1(1,1)^2+k1(1,2)^2+k1(2,1)^2+k1(2,2)^2)/2;
    if size_kchange<taok%taok为迭代终止最小正阈值
        fprintf("Termination condition |kasi*deltak| = %d < taok =%d met. End\n", size_kchange, taok);
        break;
    else
        K=K+kasi*deltak;
    end
    
    disp(['  Optimal K found at i=', num2str(i), '  episodes is']);
    disp(K);
    % disp(' replace K with K01');
    % K=K01;
    % disp(K);
end  % end for K
 
 disp(['Optimal K is found after ', num2str(i), '  episodes']);
 format long; disp(K); format short;
 disp(" ****  Done. Opimize K for Mimimum Variance  **** ");   
 
 % check the stability
[stableFlag,eigAnetc,numEigOne] = chkEigAc(A,B,K,L);
if (stableFlag==false)
    warning("Optimization gives an Un-stable control gain matrix");
end

   nID=input('Would you like to view with the eigenvalues? \n','s');
    if nID=='y' ||  nID=='Y'      
        DEBUGPLOT=true;
    else
        return;
    end
 
 if DEBUGPLOT==true
     figure('name','K');
     plot(allK(:,1,1),'.-b'); hold on;
     plot(allK(:,1,2),'.-r');
     plot(allK(:,2,1),'o-b');
     plot(allK(:,2,2),'o-r');
     title('gain matrix K (blue: k11,k21, red k12,k22)');
     xlabel('training episode');
     
     figure('name','Ycopt');
     plot(J1x,'-ob','LIneWidth',2); title('J1=Ycopt'); xlabel('main iteration i');
     % pause(0.2);
     
     chkEigAc(A,B,K,L,1);
 end

% [tmin, feas] = feasp(lmis,options)    