%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This LMI is for the paper 'Chang2020' (Theorem 1), it can be 
%              used to solve the problem of dynamic controller             
%
%              This function is based on LMI_Stability_Chang2014_case_2.m
% Author:      Yan Zong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [gamma, K] = Optim_fun(alpha, beta)
% function Gamma = Optim_fun(theta, epsilon)
    
    theta=7.9201
    epsilon=0.1499

    %%% system parameters
    A=[1 1; 0 1];
    B=[1 0; 0 1];
    E=[1 0 0 0 -1; 0 1 0 0 0];

    C1=[1 0];
    F=[0 0 1 0 0];
    D=zeros(1,2);

    C2=[1 0;0 1];
    H=[0 0 1 0 0;0 0 0 1 0];

    n_x=2;    % x的维数为2*1
    m_u=2;    % u的维数为2*1
    q_z=1;    % z的维数为1*1
    f_y=2;    % y的维数为2*1
    v_w=5;    % w的维数为5*1

    % adjusting matrix N^(n+v,n+q): three cases
    % Case 1: if v>q, N=[eye(n_x+q_z); zeros(v_w-q_z,n_x+q_z)]
    N1=eye(2); % 2x2 
    N2=zeros(2); % 2x2
    N3=zeros(2,1); % 2x1
    N4=zeros(2); % 2x2
    N5=eye(2); % 2x2
    N6=[0;0]; % 2x1
    N7=[zeros(1,2)]; % 1x2
    N8=[zeros(1,2)]; % 1x2
    N9=eye(1); % 1x1
    N=[N1 N2 N3; N4 N5 N6; N7 N8 N9];

    %-------Initial a LMI system--------
    setlmis([]);   
    % U=[U1 U2; U3 U4]
    U1=lmivar(2,[2 2]); % 2x2 matrix 
    U2=lmivar(2,[2 2]); % 2x2 matrix 
    U3=lmivar(2,[2 2]); % 2x2 matrix 
    U4=lmivar(2,[2 2]); % 2x2 matrix 

    % V=[V1 V2; V3 V4]
    V1=lmivar(2,[2 2]); % 2x2 matrix 
    V2=lmivar(2,[2 2]); % 2x2 matrix 
    V3=lmivar(2,[2 2]); % 2x2 matrix 
    V4=lmivar(2,[2 2]); % 2x2 matrix    
    
    % P=[P1 P2; P3 P4]
    P1=lmivar(2,[2 2]); % 2x2 matrix 
    P2=lmivar(2,[2 2]); % 2x2 matrix 
    P3=lmivar(2,[2 2]); % 2x2 matrix 
    P4=lmivar(2,[2 2]); % 2x2 matrix         
    
    % G=[G1 G2 G3; G4 G5 G6; G7 G8 G9]
    G1=lmivar(2,[2 2]); % 2x2 matrix 
    G2=lmivar(2,[2 2]); % 2x2 matrix
    G3=lmivar(2,[2 1]); % 1x1 matrix
    G4=lmivar(2,[2 2]); % 2x2 matrix
    G5=lmivar(2,[2 2]); % 2x2 matrix 
    G6=lmivar(2,[2 1]); % 1x1 matrix
    G7=lmivar(2,[1 2]); % 1x1 matrix
    G8=lmivar(2,[1 2]); % 1x1 matrix
    G9=lmivar(2,[1 1]); % 1x1 matrix    
    
    gamma=lmivar(1,[1 1]);  % 1x1 symmetric matrix % performance index

    %----------L=1 1*1---------
    % the first LMI 
    lmiterm([1 1 1 P1],-1,1);   
    lmiterm([1 1 1 P2],-1,1);   
    lmiterm([1 1 1 P3],-1,1);   
    lmiterm([1 1 1 P4],-1,1);                                      

    %----------L=2 8*8---------
    % the second LMI
    % (1, 1)
    lmiterm([2 1 1 P1],-1,1);
    lmiterm([2 1 1 V4],epsilon*N1*B,C2,'s');    
    lmiterm([2 1 1 V2],epsilon*N2,C2,'s');
    lmiterm([2 1 1 V4],epsilon*N3*D,C2,'s');
    
    % (2, 1)
    lmiterm([2 2 1 P3],-1,1);        
    lmiterm([2 2 1 V4],epsilon*N4*B,C2);
    lmiterm([2 2 1 V2],epsilon*N5,C2);
    lmiterm([2 2 1 V4],epsilon*N6*D,C2);    
    lmiterm([2 1 2 V3],epsilon*N1*B,1);
    lmiterm([2 1 2 V1],epsilon*N2,1);
    lmiterm([2 1 2 V3],epsilon*N3*D,1);    
    
    % (2, 2)
    lmiterm([2 2 2 P4],-1,1);
    lmiterm([2 2 2 V3],epsilon*N4*B,1,'s');
    lmiterm([2 2 2 V1],epsilon*N5,1,'s');
    lmiterm([2 2 2 V3],epsilon*N6*D,1,'s');
    
    % (3, 1)    
    lmiterm([2 3 1 V4],epsilon*N7*B,C2);
    lmiterm([2 3 1 V2],epsilon*N8,C2);    
    lmiterm([2 3 1 V4],epsilon*N9*D,C2);    
    lmiterm([2 1 3 V4],epsilon*N1*B,H);
    lmiterm([2 1 3 V2],epsilon*N2,H);    
    lmiterm([2 1 3 V4],epsilon*N3*D,H);    
            
    % (3, 2)
    lmiterm([2 3 2 V3],epsilon*N7*B,1);
    lmiterm([2 3 2 V1],epsilon*N8,1);    
    lmiterm([2 3 2 V3],epsilon*N9*D,1);    
    lmiterm([2 2 3 V4],epsilon*N4*B,H);
    lmiterm([2 2 3 V2],epsilon*N5,H);    
    lmiterm([2 2 3 V4],epsilon*N6*D,H);   
    
    % (3, 3)
    lmiterm([2 3 3 gamma],-1,1);        
    lmiterm([2 3 3 V4],epsilon*N7*B,H,'s');
    lmiterm([2 3 3 V2],epsilon*N8,H,'s');
    lmiterm([2 3 3 V4],epsilon*N9*D,H,'s');    
    
    % (4, 1)
    lmiterm([2 4 1 G1],1,A);
    lmiterm([2 4 1 G3],1,C1);
    lmiterm([2 4 1 V4],B,C2);
    
    % (4, 2)
    lmiterm([2 4 2 V3],B,1);
    
    % (4, 3)
    lmiterm([2 4 3 G1],1,E);
    lmiterm([2 4 3 G3],1,F);
    lmiterm([2 4 3 V4],B,H);
    
    % (4, 4)
    lmiterm([2 4 4 G1],-1,1,'s');        
    lmiterm([2 4 4 P1],1,1);
    
    % (5, 1)
    lmiterm([2 5 1 G4],1,A);
    lmiterm([2 5 1 G6],1,C1);
    lmiterm([2 5 1 V2],1,C2);
    
    % (5, 2)
    lmiterm([2 5 2 V1],1,1);
    
    % (5, 3)
    lmiterm([2 5 3 G4],1,E);
    lmiterm([2 5 3 G6],1,F);
    lmiterm([2 5 3 V2],1,H);
    
    % (5, 4)
    lmiterm([2 5 4 G4],-1,1);
    lmiterm([2 4 5 G2],-1,1);
    lmiterm([2 5 4 P3],1,1);
    
    % (5, 5)
    lmiterm([2 5 5 G5],-1,1,'s');
    lmiterm([2 5 5 P4],1,1);    
    
    % (6, 1)
    lmiterm([2 6 1 G7],1,A);
    lmiterm([2 6 1 G9],1,C1);
    lmiterm([2 6 1 V4],D,C2);
    
    % (6, 2)
    lmiterm([2 6 2 V3],D,1);
    
    % (6, 3)
    lmiterm([2 6 3 G7],1,E);
    lmiterm([2 6 3 G9],1,F);
    lmiterm([2 6 3 V4],D,H);
    
    % (6, 4)
    lmiterm([2 6 4 G7],-1,1);
    lmiterm([2 4 6 G3],-1,1);
    
    % (6, 5)
    lmiterm([2 6 5 G8],-1,1);
    lmiterm([2 5 6 G6],-1,1);
    
    % (6, 6)
    lmiterm([2 6 6 G9],-1,1,'s');
    lmiterm([2 6 6 0],1);
    
    % (7, 1)
    lmiterm([2 7 1 V2],1,C2);
    lmiterm([2 1 7 U3],-epsilon*N1*B,1);
    lmiterm([2 1 7 U1],-epsilon*N2,1);
    lmiterm([2 1 7 U3],-epsilon*N3*D,1);
    
    % (7, 2)
    lmiterm([2 7 2 V1],1,1);
    lmiterm([2 2 7 U3],-epsilon*N4*B,1);
    lmiterm([2 2 7 U1],-epsilon*N5,1);
    lmiterm([2 2 7 U3],-epsilon*N6*D,1);
    
    % (7, 3)
    lmiterm([2 7 3 V2],1,H);
    lmiterm([2 3 7 U3],-epsilon*N7*B,1);
    lmiterm([2 3 7 U1],-epsilon*N8,1);
    lmiterm([2 3 7 U3],-epsilon*N9*D,1);
    
    % (7, 4)
    lmiterm([2 4 7 G2],theta,1);
    lmiterm([2 4 7 U3],B,-1);
    
    % (7, 5)
    lmiterm([2 5 7 G5],theta,1);
    lmiterm([2 5 7 U1],-1,1);
    
    % (7, 6)
    lmiterm([2 6 7 G8],theta,1);
    lmiterm([2 6 7 U3],D,-1);
    
    % (7, 7)
    lmiterm([2 7 7 U1],-1,1,'s');
    
    % (8, 1)
    lmiterm([2 8 1 V4],(B'*B+D'*D),C2);    
    lmiterm([2 1 8 U4],-epsilon*N1*B,1);
    lmiterm([2 1 8 U2],-epsilon*N2,1);
    lmiterm([2 1 8 U4],-epsilon*N3*D,1);
        
    % (8, 2)
    lmiterm([2 8 2 V3],(B'*B+D'*D),1);    
    lmiterm([2 2 8 U4],-epsilon*N4*B,1);
    lmiterm([2 2 8 U2],-epsilon*N5,1);
    lmiterm([2 2 8 U4],-epsilon*N6*D,1);
    
    % (8, 3)
    lmiterm([2 8 3 V4],(B'*B+D'*D),H);    
    lmiterm([2 3 8 U4],-epsilon*N7*B,1);
    lmiterm([2 3 8 U2],-epsilon*N8,1);
    lmiterm([2 3 8 U4],-epsilon*N9*D,1);    
    
    % (8, 4)
    lmiterm([2 4 8 G1],theta,B);
    lmiterm([2 4 8 G3],theta,D);
    lmiterm([2 4 8 U4],B,-1);
        
    % (8, 5)
    lmiterm([2 5 8 G4],theta,B);
    lmiterm([2 5 8 G6],theta,D);
    lmiterm([2 5 8 U2],-1,1);
    
    % (8, 6)
    lmiterm([2 6 8 G7],theta,B);
    lmiterm([2 6 8 G9],theta,D);
    lmiterm([2 6 8 U4],D,-1);    
    
    % (8, 7)
    lmiterm([2 8 7 U3],(B'*B+D'*D),-1); 
    lmiterm([2 7 8 U2],1,-1); 
    
    % (8, 8)
    lmiterm([2 8 8 U4],(B'*B+D'*D),-1,'s');
    
    %%%%%%%%%%%优化gamma%%%%%%%%%%%%%%%%
    lmis=getlmis;
    n=decnbr(lmis); % the number of decision variables in LMI system
    c=zeros(n,1); % n x 1 zeros
    for j=1:n
    [gaj]=defcx(lmis,j,gamma); % returns [gaj] of the matrix variables with gamma 
                               % when the j-th decision variable is set to 1 
                               % and all others to 0

    c(j)=trace(gaj);  % sum of diagonal elements.
    end
    options=[1e-5,500,0,0,0];
    [copt,xopt]=mincx(lmis,c,options);

    Gamma=sqrt(dec2mat(lmis,xopt,gamma));
    
    u1=dec2mat(lmis,xopt,U1);
    u2=dec2mat(lmis,xopt,U2);
    u3=dec2mat(lmis,xopt,U3);
    u4=dec2mat(lmis,xopt,U4);
    u=[u1 u2; u3 u4];
    
    v1=dec2mat(lmis,xopt,V1);
    v2=dec2mat(lmis,xopt,V2);
    v3=dec2mat(lmis,xopt,V3);
    v4=dec2mat(lmis,xopt,V4);    
    v=[v1 v2; v3 v4];
    
    K=theta*inv(u)*v;    

%     clear lmis, su1, su2, u1, u2, U;
%     clear lmis, sv1, sv2, v1, v2, V;
%     clear G1, G2, G3, G4, P, gamma;
    
%     Gamma    
%     K
    
% end