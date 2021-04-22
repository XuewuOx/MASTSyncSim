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
% clear all;
% close all;
% clc;


%% Simulaiton Configuration 1: Network Topology
disp("Clock Synchronization Simulation");

szsim=800;   %总时间
T=1; %clock synchronization interval
arbitraryNetwork=false; % false; %true;


if arbitraryNetwork 
    % Random network generation
    nNode=20; nEdge=20;  
%     nNode=4; nEdge=4;  
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
        ];%单线
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
        ];%环形
    % Lsu12=[2 -1 -1 0 0 0 0 0 0 0 0 0;
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
    % ];%树形
    L=L12c;
    % L=L4;
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

% (1) genertion of process and measurement noises
% According to the 3Sigma rule of Gaussian distribution
%  99.73% of the noise will be in the range of [mu-3*sigma, mu+3*sigma]
% This suggests nearly all values are taken to lie within 
% three standard deviations of the mean, and thus it is 
% empirically useful to treat 99.7% probability
%  as near certainty.
% Therefore, the noise range is -/+ 3*sigma
sigma1sqr=10^-12;%有噪声 variance of offset noise 
sigma2sqr=10^-12;  % variance of skew noise
sigma3sqr=10^-12; % variance of offset's measurement noise, 
               % this is equal to the noise caused by the delay jitters
% TODO 噪声的方差 sigma^2 可根据WSN硬件平台实测数据进行调整
%  根据Zong的实测数据， CPU处理延时的均值是311us,
%                       CPU处理延时的抖动(standard deviation)是4us
%  CPU处理延时相当于在执行器上的扰动，理论上需要根据控制器增益K，
%  通过K^-1 换算为 系统输出的测量值扰动  
%  这里，暂且将CPU处理延迟的抖动作为测量噪声， 因此 sigma3=4us
 sigma3sqr=(4*10^(-6))^2;  % 1.6*10^-11
               
% sigma=sigma3sqr*[1 1;1 2];
mu=[0 0];
R=[sigma1sqr 0;0 sigma2sqr];%过程噪声协方差矩阵
% Q1=[Q zeros(2,2);zeros(2,2) Q];
Q=[sigma3sqr sigma3sqr;sigma3sqr 2*sigma3sqr];%测量噪声协方差矩阵
%  R1=[R zeros(2,2);zeros(2,2) R]
fprintf("Noise levels:\r");
fprintf("    Phase (offset) noise std = %d\r", sqrt(sigma1sqr));
fprintf("    Freq (skew) noise std = %d\r", sqrt(sigma2sqr));
fprintf("    Observation (offset) noise std = %d\r", sqrt(sigma3sqr));

% Q=eye(2)

% for k-th node, the noises is a 4-by-sz matrix stated in  
% odd row of procNoise are (1) process noses of theta  (1-by-sz) 过程噪声中thet的随机值
% even row of procNoise are (2) process noise of gamma (1-by-sz)  过程噪声中gamma的随机值
% odd row of measNoise are  (3) measurement noise of theta (1-by-sz) 测量噪声中theta 的随机值(生成多维正态数据）
% odd row of measNoise are (4) measurement noise of gamma (1-by-sz)
% 测量噪声中gamma的随机值(生成多维正态数据）
% DXW: gamma的测量噪声来自于theta的测量噪声，不应该是独立的 ????
%
% for循环实现噪声向量随节点个数自动生成
procNoisew=[];
measNoisev=[];
for k=1:nNode
procNoisew((k-1)*2+1:(k-1)*2+2,:)=[sqrt(sigma1sqr)*randn(1,szsim);sqrt(sigma2sqr)*randn(1,szsim)];
measNoisev((k-1)*2+1:(k-1)*2+2,:)=mvnrnd(mu,Q,szsim)'; 
end


% (3) TODO: Controller design
% D =[0.377684488137080   0.046252370994938;0.030384079824661   0.384554203344394];%线性搜索LMI 改变A值所得的K
% D=[0.4563 0.0035;0.0112 0.4723];%0.1传统LMI所得K值
% K=[0.3200 0.0450;0.2330 0.0511];%线性搜索LMI所得K值, Hu,Sec 5, 

% gain matrix for tree20 network 
   K01=[0.1 0;0 0.1];%原来的初始值
    K02=[0.3517 -0.1997;0.0159 0.2565];%以传统中得到的优化后的增益作为初始值,J变大了10-14
    K03=[0.3338 -0.0840;0.0133 0.3034];%以传统中得到的优化后的增益作为初始值 10-14  10-16
    
% K=K03;   
% K0=K;
% Linearsearch20shu
% load goodK_2069e_10.mat
% load goodK_9253e_10.mat
% load goodK_nonLMI_2698e_08.mat
% szsim=400; 
 
%  K=[0.3338 -0.0840;0.0133 0.3034];
%     mean value of all nodes at end: [offsetba=8.398589e-02 skewbar=2.098934e-04]
%    with standard deviation [3.662655e-07, 1.559589e-07]

%  K=[-0.0116 0.0014;0.0101 0.1002];%绾挎?ф悳绱㈢殑
%   mean value of all nodes at end: [offsetba=8.399398e-02 skewbar=2.100283e-04]
%    with standard deviation [2.623130e-04, 2.732936e-06]
% K=[0.0446 0.0337;0.0327 0.0562];
%   mean value of all nodes at end: [offsetba=8.398976e-02 skewbar=2.100030e-04]
%    with standard deviation [3.793082e-05, 6.394281e-06]
% K=[0.3322 -0.0797;0.0118 0.3035];
%    mean value of all nodes at end: [offsetba=8.398261e-02 skewbar=2.099536e-04]
%    with standard deviation [3.613055e-07, 1.946393e-07]
format long   
fprintf("Control Gain K="); disp(K);
format short

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



%
%  [alpha,beta]=meshgrid(0.01:0.01:0.58,0.01:0.01:0.58);
%inv(B)
%     D=[alpha 0;0 beta];
%   D=inv(B)*A
%  eig(D)

% internal variables for clock state, output and sync errors
% %noKalman
y1=zeros(2*nNode,szsim); % output, a matrix for outputs at all sz simulation steps
xx1=zeros(2*nNode,szsim); % state, a matrix for states at all sz simulation steps
y2=zeros(2*nNode,szsim); % output
xx2=zeros(2*nNode,szsim); % state

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
% disp(L);
listnonServoClk=[];
listRefClk=[];
while true
    nID=input('Any clock node you want to change to non-sevo clock? \n Type the node ID (0 for nothing to change): ');
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
yerr(:,1)=yk-kron(ones(nNode,1),Ybar(:,1));%误差向量（输出-输出均值）胡equ29  ??
   
%%  Simulating the Networked Synchronization controller 
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


A1=kron(eye(nNode),A);%克罗内克积
% D=[alpha(i,j) 0;0 beta(i,j)];
BK=B*K;  % variable D for the feedback gain matrix K in Hu2019
BK1=kron(L,BK);

% D=[alpha(i,j) 0 0 0;0 beta(i,j) 0 0;0 0 alpha(i,j) 0;0 0 0 beta(i,j)];
%    Q=[theta1(k) 0;0 theta2(k)];R=[theta3(k) theta3(k);theta3(k) theta4(k)];
    
B1=kron(eye(nNode),B*K);
L1=kron(L,eye(2));

if chkEigAc(A,B,K,L)==false
    % the close form of the NCS has unstable eigenvalue
    warning("Unstable eigenvalue of the networked closed loop system \n");
else
    disp("stable system (two eigenvalue ==1, others <1)");
end
% Please note, B1*L1 shoulde be the same as BK1
% sz=500;
for k = 2:szsim

%   % x-based simulation, State and output updates, Hu2019, eq.26, 27
    U=L1*y2(:,k-1); % get the output differece with neighbours
    % U=L1*zeros(2*nNode,1); % free running clock, no regulation
    xx2(:,k)=A1*xx2(:,k-1)-B1*U+procNoisew(:,k-1); % state updates
    y2(:,k)=xx2(:,k)+measNoisev(:,k); % output updates

    % y-base sys, noise-free      
    y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1);%w输出Y值 胡 eq 28
    xx1(:,k)=y1(:,k);%状态X值

%    % noisy y-base system
%     v1=measNoisev(:,k-1);  % DXW: replace he next two lines
%     v2=measNoisev(:,k);
%               
%     % Rearrange the closed system as Y being the state variable of interest
%    % y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1)+procNoisew(:,k-1)+v2-A1*v1;%w输出Y值 胡 eq 28
%     y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1)+procNoisew(:,k-1)+v2-A1*v1+BK1*v1;%w输出Y值 胡 eq 28
%   
%     xx1(:,k)=y(:,k)-v2;%状态X值
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
    yerr(:,k)=yk-kron(ones(nNode,1),Ybar(:,k));%误差向量（输出-输出均值）胡equ29  ??
   
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
if (length(listRefClk)==0)
   title(sprintf('Topology: %d nodes, %d edges (no refence clk)',nNode, nEdge)); 
else
   title( "Topology: "+nNode+" nodes  "+nEdge+ "edges  "+ num2str(length(listRefClk)+ " Refence clock = [" + num2str(listRefClk))+"]");
end 
subplot(1,2,2)
plot(1:k,Ystd, 'MarkerSize', 5, 'LineWidth', 2)
title('std of offset and skew among nodes');
xlabel('time (s)');
legend('std of offset among nodes','std of skew among nodes');
grid on;
fprintf("Simulation Done!  %d nodes, %d sync steps\n",nNode,szsim);
 disp("     " + num2str(length(listnonServoClk)+ " non-servo clock = [" + num2str(listnonServoClk))+"]");
 disp("     " + num2str(length(listRefClk)+ " Refence clock = [" + num2str(listRefClk))+"]");
fprintf("   mean value of all nodes at end: [offsetba=%d skewbar=%d]\n",Ybar(1,szsim), Ybar(2,szsim));
fprintf("   with standard deviation [%d, %d]\n", Ystd(1,k), Ystd(2,k));

xvar=var(xx2');
xthetavar=[1:nNode; xvar(1:2:2*nNode-1)];
xskewvar=[1:nNode; xvar(2:2:2*nNode)];

xvarsorted=sortrows(xthetavar',[2]); % sort in ascending order based on the elements in the 2 column.
disp("node (in ascending theta var) ="); disp(xvarsorted(:,1)');
disp("                    theta var ="); disp(xvarsorted(:,2)');

return;
%%  plot results for visualization
figure('name','Sychronization offset error');
hold on; cla; plot(1:k,yerr(indTheta,1:k), '-.');


y=y2;
xx=xx2;
%         figure(1)
%         plot([1:sz],xx)
%         legend('节点1 offset状态值','节点1 skew状态值','节点2 offset状态值','节点2 skew状态值','节点3 offset状态值','节点3 skew状态值','节点4 offset状态值','节点4 skew状态值');       
%         xlabel('n');ylabel('状态值');
%         
%         figure(2)
%         plot([1:sz],y)
%         legend('节点1 offset输出值','节点1 skew输出值','节点2 offset输出值','节点2 skew输出值','节点3 offset输出值','节点3 skew输出值','节点4 offset输出值','节点4 skew输出值');
%         xlabel('n');ylabel('输出值');
        
        figure(3)
        subplot(2,2,1)
        for kk=1:2:2*nNode-1
            plot([1:szsim],xx(kk,:));hold on
        end
      %    legend('节点1 offset状态值','节点2 offset状态值','节点3 offset状态值','节点4 offset状态值','节点5 offset状态值','节点6 offset状态值','节点7 offset状态值','节点8 offset状态值','节点9 offset状态值','节点10 offset状态值','节点11 offset状态值','节点12 offset状态值');
        xlabel('n');ylabel('状态值');
        
        subplot(2,2,2)
        for kk=2:2:2*nNode
            plot([1:szsim],xx(kk,:));hold on
        end
       % legend('节点1 skew状态值','节点2 skew状态值','节点3 skew状态值','节点4 skew状态值','节点5 skew状态值','节点6 skew状态值','节点7 skew状态值','节点8 skew状态值','节点9 skew状态值','节点10 skew状态值','节点11 skew状态值','节点12 skew状态值');
        xlabel('n');ylabel('状态值');
        
        subplot(2,2,3)
      for kk=1:2:2*nNode-1
            plot([1:szsim],y(kk,:));hold on
        end
       % legend('节点1 offset输出值','节点2 offset输出值','节点3 offset输出值','节点4 offset输出值','节点5 offset输出值','节点6 offset输出值','节点7 offset输出值','节点8 offset输出值','节点9 offset输出值','节点10 offset输出值','节点11 offset输出值','节点12 offset输出值');
        xlabel('n');ylabel('输出值');
        
        subplot(2,2,4)
        for kk=2:2:2*nNode
            plot([1:szsim],y(kk,:));hold on
        end
        % legend('节点1 skew输出值','节点2 skew输出值','节点3 skew输出值','节点4 skew输出值','节点5 skew输出值','节点6 skew输出值','节点7 skew输出值','节点8 skew输出值','节点9 skew输出值','节点10 skew输出值','节点11 skew输出值','节点12 skew输出值');
        xlabel('n');ylabel('输出值');
      
      % 
        for kk=1:2*nNode
            var_xx(kk)=var(xx(kk,:));
        end
      mean_xx=mean(xx,2);
      

        
        


        
   
