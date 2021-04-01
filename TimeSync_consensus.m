% % ***************************************************************************
% % * 
% * Created on:  31 03 2021
% * Modified on: 31.03.2021
% * 
% * Author:      Xuweu Dai
% *
% * File:        TimeSync_consensus.m
% *
% *
% * This program is free software: you can redistribute it and/or modify
% * it under the terms of the GNU Lesser General Public License as published by
% * the Free Software Foundation, either version 3 of the License, or
% * (at your option) any later version.
% 
% * This program is distributed in the hope that it will be useful,
% * but WITHOUT ANY WARRANTY; without even the implied warranty of
% * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% * GNU Lesser General Public License for more details.
% *
% %              _
% %  \/\ /\ /   /  * _  '
% % _/\ \/\/ __/__.'(_|_|_
% % **************************************************************************/
%%
clear all;
% close all;
% clc;


%% Simulaiton Configuration 1: Network Topology
disp("Clock Synchronization Simulation");

sz=200;   %��ʱ��
T=1; %clock synchronization interval
arbitraryNetwork=true; % false; %true;


if arbitraryNetwork 
    % Random network generation
    nNode=20; nEdge=20;  
    [G0 Tree0]=genNet(nNode,nEdge,1);
    if isempty(G0)
       disp('Network generation and simlation cancelled');
       return;
    end
    disp('Do you want to use the graph or spanning tree for TimeSync?');
    disp('     g/G -- use graph;   t/T --- use the minimum spannning tree');
    usrinput = input( 'Type your choice: ', 's');
    if usrinput == 't' || usrinput=='T' || usrinput=='n' || usrinput=='N'
        fprintf("the TREE is used for time sync. %d nodes, %d edges\n", ...
            nNode,numedges(Tree0));
        netG=Tree0;
    else
        fprintf("the GRAPH is used for time sync. %d nodes, %d edges\n", ...
            nNode,numedges(G0));
        netG=G0;
    end
    L=full(laplacian(netG));
    
else
    % Network topology specified by the Laplacian matrix L
    L3=[1 -1 0;
        -1 2 -1;
        0 -1 1];
    L4=[1 -1 0 0;
        -1 2 -1 0;
        0 -1 2 -1;
        0 0  -1 1];
    
    L4r=[1 -1 0 0;
        -1 3 -1 -1;
        0 -1 1 0;
        0 -1 0 1];
    % L12: 12 nodes linear network
    L12=[1 -1 0 0 0 0 0 0 0 0 0 0;
        -1 2 -1 0 0 0 0 0 0 0 0 0;
        0 -1 2 -1 0 0 0 0 0 0 0 0;
        0 0 -1 2 -1 0 0 0 0 0 0 0;
        0 0 0 -1 2 -1 0 0 0 0 0 0;
        0 0 0 0 -1 2 -1 0 0 0 0 0;
        0 0 0 0 0 -1 2 -1 0 0 0 0;
        0 0 0 0 0 0 -1 2 -1 0 0 0;
        0 0 0 0 0 0 0 -1 2 -1 0 0;
        0 0 0 0 0 0 0 0 -1 2 -1 0;
        0 0 0 0 0 0 0 0 0 -1 2 -1;
        0 0 0 0 0 0 0 0 0 0 -1 1;
        ];%����
    % L12c: 12 nodes circular network
    L12c=[2 -1 0 0 0 0 0 0 0 0 0 -1;
        -1 2 -1 0 0 0 0 0 0 0 0 0;
        0 -1 2 -1 0 0 0 0 0 0 0 0;
        0 0 -1 2 -1 0 0 0 0 0 0 0;
        0 0 0 -1 2 -1 0 0 0 0 0 0;
        0 0 0 0 -1 2 -1 0 0 0 0 0;
        0 0 0 0 0 -1 2 -1 0 0 0 0;
        0 0 0 0 0 0 -1 2 -1 0 0 0;
        0 0 0 0 0 0 0 -1 2 -1 0 0;
        0 0 0 0 0 0 0 0 -1 2 -1 0;
        0 0 0 0 0 0 0 0 0 -1 2 -1;
        -1 0 0 0 0 0 0 0 0 0 -1 2;
        ];%����
    % L=[2 -1 -1 0 0 0 0 0 0 0 0 0;
    %   -1 3 0 -1 -1 0 0 0 0 0 0 0;
    %   -1 0 3 0 0 -1 -1 0 0 0 0 0;
    %   0 -1 0 3 0 0 0 -1 -1 0 0 0;
    %   0 -1 0 0 3 0 0 0 0 -1 -1 0;
    %   0 0 -1 0 0 1 0 0 0 0 0 0;
    %   0 0 -1 0 0 0 2 0 0 0 0 -1;
    %   0 0 0 -1 0 0 0 1 0 0 0 0;
    %   0 0 0 -1 0 0 0 0 1 0 0 0;
    %   0 0 0 0 -1 0 0 0 0 1 0 0;
    %   0 0 0 0 -1 0 0 0 0 0 1 0;
    %   0 0 0 0 0 0 -1 0 0 0 0 1
    % ];%����
    L=L12;
    [netG,L]=genNetbyL();
    % nNode=length(L);% number of nodes
    % nEdge=trace(L)/2; % number of edges
end

% update the nNode and nEdge, as they may be changed due to spannning tree
nNode=numnodes(netG);% number of nodes
nEdge=numedges(netG); % number of edges
% any changes to the network and clock configuraiton 
disp(sprintf("Network created with : %d nodes.%d edges", nNode,nEdge));
% disp(L);
%% Simulaiton Configuration 2: Clock, noises and controller

A=[1 T;0 1];B=[1 T;0 1];H=[1 0;0 1];
B2=zeros(4,4);
% D =[0.377684488137080   0.046252370994938;0.030384079824661   0.384554203344394];%��������LMI �ı�Aֵ���õ�K
% D=0;
% D=[0.3200 0.0450;0.2330 0.0511];%��������LMI����Kֵ, Hu,Sec 5, 
% D=[0.4563 0.0035;0.0112 0.4723];%0.1��ͳLMI����Kֵ
K=[0.3200 0.0450;0.2330 0.0511];%��������LMI����Kֵ, Hu,Sec 5, 
fprintf("Control Gain K="); disp(K);


% (1) Initial value of Clock offset and skew
% (1a) set initial state values
% 1.  Identical clocks: all nodes's [offset, skew]=[100, 20]
x0_noisefree=10^(-6)*kron(ones(nNode,1),[100,20]');
% 10^(-6)*[100; 20; 100; 20; 100; 20; 100; 20; 100; 20; 100; 20; 100; 20; ...
%            100; 20; 100; 20; 100; 20; 100; 20; 100; 20];

% 2. Non-identical clocks
a0=500; % offset 100us 
b0=20; % skew 20ppm
a=40; % offset inittial value variation range: +/- 40 uniform distribution
b=20; % skew initial value variation raange: +/- 20 uniform distribution
iniOffset=a0+a*2*(rand(1,nNode)-0.5);
iniSkew=b0+b*2*(rand(1,nNode)-0.5);
x0=10^(-6)*reshape([iniOffset;iniSkew],[],1);
x0=10^(-6)*[10:10:2*nNode*10];
% x0(2)=0; x0(1)=0; % set node 1 as ideal pefect clock
% x0=10^(-6)*[99.36; 19.25; 99.3; 19.3; 98.7; 18.5; 99.8; 19.8; 99.36; 19.5; ...
%             99.3; 19.3; 98.7; 18.5; 99.8; 19.8; 99.36;19.5; 99.3;  19.3;
%             98.7; 18.5; 99.8; 19.8];
fprintf("Clock initial values:\r");
fprintf("    Offset %d +/- %d us \r", a0, a);
fprintf("    Skew  %d +/- %d ppm\r ", b0, b);



% (2) genertion of process and measurement noises
% According to the 3Sigma rule of Gaussian distribution
%  99.73% of the noise will be in the range of [mu-3*sigma, mu+3*sigma]
% This suggests nearly all values are taken to lie within 
% three standard deviations of the mean, and thus it is 
% empirically useful to treat 99.7% probability
%  as near certainty.
% Therefore, the noise range is -/+ 3*sigma
sigma1sqr=10^-12;%������ variance of offset noise 
sigma2sqr=10^-12;  % variance of skew noise
sigma3sqr=10^-12; % variance of offset's measurement noise, 
               % this is equal to the noise caused by the delay jitters
% TODO �����ķ��� sigma^2 �ɸ���WSNӲ��ƽ̨ʵ�����ݽ��е���
%  ����Zong��ʵ�����ݣ� CPU������ʱ�ľ�ֵ��311us,
%                       CPU������ʱ�Ķ���(standard deviation)��4us
%  CPU������ʱ�൱����ִ�����ϵ��Ŷ�����������Ҫ���ݿ���������K��
%  ͨ��K^-1 ����Ϊ ϵͳ����Ĳ���ֵ�Ŷ�  
%  ������ҽ�CPU�����ӳٵĶ�����Ϊ���������� ��� sigma3=4us
 sigma3sqr=(4*10^(-6))^2*0;  % 1.6*10^-11
               
sigma=sigma3sqr*[1 1;1 2];
mu=[0 0];
Q=[sigma1sqr 0;0 sigma2sqr];%��������Э�������
% Q1=[Q zeros(2,2);zeros(2,2) Q];
R=[sigma3sqr sigma3sqr;sigma3sqr 2*sigma3sqr];%��������Э�������
%  R1=[R zeros(2,2);zeros(2,2) R]
fprintf("Noise levels:\r");
fprintf("    Phase (offset) noise std = %d\r", sqrt(sigma1sqr));
fprintf("    Freq (skew) noise std = %d\r", sqrt(sigma2sqr));
fprintf("    Obsrvation (offset) noise std = %d\r", sqrt(sigma3sqr));

% Q=eye(2)

% for k-th node, the noises is a 4-by-sz matrix stated in  
% odd row of procNoise are (1) process noses of theta  (1-by-sz) ����������thet�����ֵ
% even row of procNoise are (2) process noise of gamma (1-by-sz)  ����������gamma�����ֵ
% odd row of measNoise are  (3) measurement noise of theta (1-by-sz) ����������theta �����ֵ(���ɶ�ά��̬���ݣ�
% odd row of measNoise are (4) measurement noise of gamma (1-by-sz)
% ����������gamma�����ֵ(���ɶ�ά��̬���ݣ�
% DXW: gamma�Ĳ�������������theta�Ĳ�����������Ӧ���Ƕ����� ????
%
% forѭ��ʵ������������ڵ�����Զ�����
procNoisew=[];
measNoisev=[];
for k=1:nNode
procNoisew((k-1)*2+1:(k-1)*2+2,:)=[sqrt(sigma1sqr)*randn(1,sz);sqrt(sigma2sqr)*randn(1,sz)];
measNoisev((k-1)*2+1:(k-1)*2+2,:)=mvnrnd(mu,sigma,sz)'; 
end

%
%  [alpha,beta]=meshgrid(0.01:0.01:0.58,0.01:0.01:0.58);
%inv(B)
%     D=[alpha 0;0 beta];
%   D=inv(B)*A
%  eig(D)

% internal variables for clock state, output and sync errors
% %noKalman
y1=zeros(2*nNode,sz); % output, a matrix for outputs at all sz simulation steps
xx1=zeros(2*nNode,sz); % state, a matrix for states at all sz simulation steps
y2=zeros(2*nNode,sz); % output
xx2=zeros(2*nNode,sz); % state

% yerr=zeros(2*nNode,sz); % output errors

indTheta=1:2:nNode*2-1;  % row indices for theta
indSkew=2:2:nNode*2; % row indices for skew

% xx1(:,1)=x0_noisefree(1:2*nNode);
xx1(:,1)=x0;
y1(:,1)=xx1(:,1)+measNoisev(:,1);
    
xx2(:,1)=x0(1:2*nNode);
y2(:,1)=xx2(:,1)+measNoisev(:,1);
xx2(:,1)=x0;
y2(:,1)=x0;

% last check before start simulation
disp(sprintf("Network size: %d nodes. Topology: L=", nNode));
disp(L);
listnonServoClk=[];
listRefClk=[];
while true
    nID=input('Any clock node you want to change to non-sevo clock? Type the node ID. 0 for nothing to change,');
    if nID==0
        break;
    end
    if (nID<0 || nID >nNode)
        disp("Invalid node ID. Please type again");
    end
    % change the nID node to non-servo clock
    %  set the nID row and col to zero vector (non-servo clock)
    %  DXW debug: should onle reset nID row to zero, nID col should be
    %  unchanged
    L(nID,:)=zeros(1,nNode); % ??? 
    % L(:,nID)=zeros(1,nNode)';
    listnonServoClk=[listnonServoClk,nID]; % add nID to the list
    fprintf(" ... ... Done. Node %d now is drifting, but non-servo clock\n",nID);

    % Check to set the nID node as reference clock?
    refClk=input(sprintf('Set node %d as a reference clock (theta=0, gamma=0, procNoisew=0) (Y) or just non-servo clock (N)? ',...
            nID),'s');
    if refClk=='y' || refClk=='Y'
        x0((nID-1)*2+1:(nID-1)*2+2)=0; % initial offset and gamme to zero
        xx2(:,1)=x0;
        y2(:,1)=x0;
        procNoisew([(nID-1)*2+1:(nID-1)*2+2],:)=0; % clear process noises w
        listRefClk=[listRefClk,nID];
            fprintf(" ... ... Done. Set node %d to reference clock.\n", nID);
    end
    fprintf("Revised Network and Clock: %d nodes.%d edges,\n", nNode,nEdge);
   disp("     " + num2str(length(listnonServoClk))+ " non-servo clock = [" + num2str(listnonServoClk)+"]");
   disp("     " + num2str(length(listRefClk))+ " Refence clock = [" + num2str(listRefClk)+"]");
   [netG L]=genNetbyL(L);
end



% initial errors for interation from k=2  
yk=y1(:,1);
yk=y2(:,1);
Ybar(:,1)=[mean(yk(indTheta)');mean(yk(indSkew)')];
yerr(:,1)=yk-kron(ones(nNode,1),Ybar(:,1));%������������-�����ֵ����equ29  ??
   
% Simulating the Networked Synchronization controller 
disp('Now start simulation ');
figure('Name','Simulation Animation'); 
        cm = colormap('Lines'); 
    subplot(2,2,1);
    title('offset \theta of all nodes & avg');
    subplot(2,2,2);
    title('skew \gamma of all nodes & avg');
    subplot(2,2,3);
    title('errors of \theta wrt.the average');
    % ax = axes; ax.ColorOrder = cm;
    subplot(2,2,4);
    title('errors of \gamma wrt.the average');
    % ax = axes; ax.ColorOrder = cm;


A1=kron(eye(nNode),A);%�����ڿ˻�
% D=[alpha(i,j) 0;0 beta(i,j)];
BK=B*K;  % variable D for the feedback gain matrix K in Hu2019
BK1=kron(L,BK);
% B2=[B1 -B1;-B1 B1]    
% D=[alpha(i,j) 0 0 0;0 beta(i,j) 0 0;0 0 alpha(i,j) 0;0 0 0 beta(i,j)];
%    Q=[theta1(k) 0;0 theta2(k)];R=[theta3(k) theta3(k);theta3(k) theta4(k)];
    
B1=kron(eye(nNode),B*K);
L1=kron(L,eye(2));
% Please note, B1*L1 shoulde be the same as BK1
% sz=500;
for k = 2:sz

%   % x-based simulation, State and output updates, Hu2019, eq.26, 27
    U=L1*y2(:,k-1); % get the output differece with neighbours
    % U=L1*zeros(2*nNode,1); % free running clock, no regulation
    xx2(:,k)=A1*xx2(:,k-1)-B1*U+procNoisew(:,k-1); % state updates
    y2(:,k)=xx2(:,k)+measNoisev(:,k); % output updates

    % y-base sys, noise-free      
    y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1);%w���Yֵ �� eq 28
    xx1(:,k)=y1(:,k);%״̬Xֵ

%    % noisy y-base system
%     v1=measNoisev(:,k-1);  % DXW: replace he next two lines
%     v2=measNoisev(:,k);
%               
%     % Rearrange the closed system as Y being the state variable of interest
%    % y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1)+procNoisew(:,k-1)+v2-A1*v1;%w���Yֵ �� eq 28
%     y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1)+procNoisew(:,k-1)+v2-A1*v1+BK1*v1;%w���Yֵ �� eq 28
%   
%     xx1(:,k)=y(:,k)-v2;%״̬Xֵ
%     % proof: how x and y are swaped
%     % y=x+v2
%     % xk+v2=A*(xk-1+v1)-BK*(xk-1+v1)+wk-1+v2-A*v1+BK*v1
%     % xk=A*(xk-1)-BK*(xk-1)+wk-1 


%     % Debug codes to check if there is any mismatch between y-based simulation and
%     % x-based one
%     if sum(abs((y2(:,k)-y1(:,k)))) ~= 0        
%         warning("Different output values at k=%d", k);
%     end
%     if sum(abs((xx2(:,k)-xx1(:,k)))) ~= 0        
%         warning("Different state values at k=%d", k);
%     end

    
    % collect the synchronization error for results analysis
    % Selet which reuslts to be used: 2 for x-based simulaiton, 1 for
    % y-based
    xx=xx2;  % x-based simulaiton results % xx=xx1; for y-based simulation 
    y=y2;  % x-based simulaiton results % Or y=y1; for y-based simulation
    
    yk=y(:,k);
    Ybar(:,k)=[mean(yk(indTheta)');mean(yk(indSkew)')];
    Ystd(:,k)=[std(yk(indTheta)');std(yk(indSkew)')];
    yerr(:,k)=yk-kron(ones(nNode,1),Ybar(:,k));%������������-�����ֵ����equ29  ??
   
    subplot(2,2,1);
    drawTrajectory(xx(indTheta,1:k))
    hold on;  plot(1:k,Ybar(1,1:k), '-k', 'MarkerSize', 5, 'LineWidth', 2);
    subplot(2,2,2);
    drawTrajectory(xx(indSkew,1:k))
    hold on;  plot(1:k,Ybar(2,1:k), '-sk', 'MarkerSize', 5,'LineWidth', 2);
    subplot(2,2,3);
    hold on; cla; plot(1:k,yerr(indTheta,1:k), '-o');
    subplot(2,2,4);
    hold on; cla; plot(1:k,yerr(indSkew,1:k), '-o','LineWidth', 1);
   
    if (k<10) || (k<100 && (mod(k,10)==0)) ||(k>100 && mod(k,50)==0)
       drawnow
    end
    %pause(0.1)

end
disp('Simulation Ends');
%% plot results
simNotes=sprintf('Clk %d, ideal ref clk',6); xpos=1;nlayer=15;
figure('Name','Network Graph and Offset variance wrt time'); 
subplot(1,2,1);
if issymmetric(L)
   plot(netG); %undirected graph
else
   % idrected graph. transpose the Laplacian for correct visulzation of
   % link's directiorn
   d=diag(L);     % D is the degree matrix 
    A=diag(diag(L))-L'; 
   plot(digraph(A));
end
text(xpos-0.2,nlayer+0.2,simNotes);
title(sprintf('Topology: %d nodes, %d edges (ref clk = [ ])',nNode, nEdge));
subplot(1,2,2)
plot(1:k,Ystd, 'MarkerSize', 5, 'LineWidth', 2)
title('std of offset and skew among nodes');
xlabel('time (s)');
legend('std of offset among nodes','std of skew among nodes');

fprintf("Simulation Done!  %d nodes, %d sync steps\n",nNode,sz);
 disp("     " + num2str(length(listnonServoClk)+ " non-servo clock = [" + num2str(listnonServoClk))+"]");
 disp("     " + num2str(length(listRefClk)+ " Refence clock = [" + num2str(listRefClk))+"]");
fprintf("   mean value of all nodes at end: [offsetba=%d skewbar=%d]\n",Ybar(1,sz), Ybar(2,sz));
fprintf("   with standard deviation [%d, %d]\n", Ystd(1,k), Ystd(2,k));

xvar=var(xx2');
xthetavar=[1:nNode; xvar(1:2:2*nNode-1)];
xskewvar=[1:nNode; xvar(2:2:2*nNode)];

xvarsorted=sortrows(xthetavar',[2]); % sort in ascending order based on the elements in the 2 column.
disp("node (in ascending theta var) ="); disp(xvarsorted(:,1)');
disp("                    theta var ="); disp(xvarsorted(:,2)');

return;
%%  plot results for visualization
y=y2;
xx=xx2;
%         figure(1)
%         plot([1:sz],xx)
%         legend('�ڵ�1 offset״ֵ̬','�ڵ�1 skew״ֵ̬','�ڵ�2 offset״ֵ̬','�ڵ�2 skew״ֵ̬','�ڵ�3 offset״ֵ̬','�ڵ�3 skew״ֵ̬','�ڵ�4 offset״ֵ̬','�ڵ�4 skew״ֵ̬');       
%         xlabel('n');ylabel('״ֵ̬');
%         
%         figure(2)
%         plot([1:sz],y)
%         legend('�ڵ�1 offset���ֵ','�ڵ�1 skew���ֵ','�ڵ�2 offset���ֵ','�ڵ�2 skew���ֵ','�ڵ�3 offset���ֵ','�ڵ�3 skew���ֵ','�ڵ�4 offset���ֵ','�ڵ�4 skew���ֵ');
%         xlabel('n');ylabel('���ֵ');
        
        figure(3)
        subplot(2,2,1)
        for kk=1:2:2*nNode-1
            plot([1:sz],xx(kk,:));hold on
        end
      %    legend('�ڵ�1 offset״ֵ̬','�ڵ�2 offset״ֵ̬','�ڵ�3 offset״ֵ̬','�ڵ�4 offset״ֵ̬','�ڵ�5 offset״ֵ̬','�ڵ�6 offset״ֵ̬','�ڵ�7 offset״ֵ̬','�ڵ�8 offset״ֵ̬','�ڵ�9 offset״ֵ̬','�ڵ�10 offset״ֵ̬','�ڵ�11 offset״ֵ̬','�ڵ�12 offset״ֵ̬');
        xlabel('n');ylabel('״ֵ̬');
        
        subplot(2,2,2)
        for kk=2:2:2*nNode
            plot([1:sz],xx(kk,:));hold on
        end
       % legend('�ڵ�1 skew״ֵ̬','�ڵ�2 skew״ֵ̬','�ڵ�3 skew״ֵ̬','�ڵ�4 skew״ֵ̬','�ڵ�5 skew״ֵ̬','�ڵ�6 skew״ֵ̬','�ڵ�7 skew״ֵ̬','�ڵ�8 skew״ֵ̬','�ڵ�9 skew״ֵ̬','�ڵ�10 skew״ֵ̬','�ڵ�11 skew״ֵ̬','�ڵ�12 skew״ֵ̬');
        xlabel('n');ylabel('״ֵ̬');
        
        subplot(2,2,3)
      for kk=1:2:2*nNode-1
            plot([1:sz],y(kk,:));hold on
        end
       % legend('�ڵ�1 offset���ֵ','�ڵ�2 offset���ֵ','�ڵ�3 offset���ֵ','�ڵ�4 offset���ֵ','�ڵ�5 offset���ֵ','�ڵ�6 offset���ֵ','�ڵ�7 offset���ֵ','�ڵ�8 offset���ֵ','�ڵ�9 offset���ֵ','�ڵ�10 offset���ֵ','�ڵ�11 offset���ֵ','�ڵ�12 offset���ֵ');
        xlabel('n');ylabel('���ֵ');
        
        subplot(2,2,4)
        for kk=2:2:2*nNode
            plot([1:sz],y(kk,:));hold on
        end
        % legend('�ڵ�1 skew���ֵ','�ڵ�2 skew���ֵ','�ڵ�3 skew���ֵ','�ڵ�4 skew���ֵ','�ڵ�5 skew���ֵ','�ڵ�6 skew���ֵ','�ڵ�7 skew���ֵ','�ڵ�8 skew���ֵ','�ڵ�9 skew���ֵ','�ڵ�10 skew���ֵ','�ڵ�11 skew���ֵ','�ڵ�12 skew���ֵ');
        xlabel('n');ylabel('���ֵ');
      
      % 
        for kk=1:2*nNode
            var_xx(kk)=var(xx(kk,:));
        end
      mean_xx=mean(xx,2);
      

        
        


        
   
