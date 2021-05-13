% % ***********************************************************************
% * Descriptoin: the code in this branch is used to simulate the network
%               under the static controller, the corresponding paper is
%               submitted to IEEE Internet of Things J. Special Issue for 
%               peer review
% * Limitation: the code can run on the MATLAB R2020b
% % ***********************************************************************
clear all;
close all;
clc;
%% Simulaiton Configuration 1: Network Topology
disp("Clock Synchronization Simulation");

szsim=800; % simulation time
arbitraryNetwork=false;
% ToDo
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
    nEdge=numedges(netG); nNode=numnodes(netG);
else
    % Network topology specified by the Laplacian matrix L
    typicalLaplacian
        
    % L=L12c;
    % L=Lhuan20; % 20 node circle network
    L=Lshu21;
    Lms=tril(L,-1)+eye(length(L)); 
    [netG,L]=genNetbyL(Lms);
   [netG,L]=genNetbyL(Lshu21);
    nNode=length(L);% number of nodes
    nEdge=trace(L)/2; % number of edges
end

% show the network topology
strnet=['Network Toplogy (', num2str(nNode), ' nodes,' num2str(nEdge), ' edges)'];
figure('name',strnet); plot(netG); title(strnet);
fprintf("Network created with : %d nodes.%d edges", nNode,nEdge);

% disp(L);
%% Simulaiton Configuration 2a: Clock and noises Q R
T=1; % clock synchronization interval
A=[1 T;0 1];
B=[1 0;0 1];
H=[1 0;0 1]; % ToDo

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

%% Simulaiton Configuration 2a: Controller Design  
% (3) TODO: Controller design
% D =[0.377684488137080   0.046252370994938;0.030384079824661   0.384554203344394];%线性搜索LMI 改变A值所得的K
% D=[0.4563 0.0035;0.0112 0.4723];%0.1传统LMI所得K值
% K=[0.3200 0.0450;0.2330 0.0511];%线性搜索LMI所得K值, Hu,Sec 5, 

Kname='K';  % give a name to the K, for saving simulation results 
                % e.g., Ka_c20, Ka_tree20
NoiseSeq=13; % name for a repoeatd simulation, for saving simulaiton reults

% 1. Controllers gain matrices for tree20 network 
   K01=[0.1 0;0 0.1];%原来的初始值
   K02=[0.3517 -0.1997;0.0159 0.2565];%以传统中得到的优化后的增益作为初始值,J变大了10-14
   K03=[0.3338 -0.0840;0.0133 0.3034];%以传统中得到的优化后的增益作为初始值 10-14  10-16

   Knewa =[0.089267697012743  -0.071498708323390;  0.097945896670515   0.299574632312155];
   Knewb=[0.136517830649692  -0.072338483097365; 0.080167639975444   0.297171905799192];
   K0=K03;
   DEBUGPLOT=false;
   Linearsearch20shu; %注释掉就是传统的
    %  K=K03;   
    %K=Knewa;
 
% 2. Controllers gain matrices for circle20 (huan20) network 
  % by using traditional LMI design (no linear search)
  % KbyLMITra_huan20
   Kc20=[0.3526 -0.1309; 0.0216  0.3234];
%   K = Kc20;
   
  
  
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




%% Simulaiton Configuration 2c: Initial states 
% (1) Initial value of Clock offset and skew
% (1a) set initial state values
% 1.  Identical clocks: all nodes's [offset, skew]=[100, 20]
x0_noisefree=10^(-6)*kron(ones(nNode,1),[100,20]');
% 10^(-6)*[100; 20; 100; 20; 100; 20; 100; 20; 100; 20; 100; 20; 100; 20; ...
%            100; 20; 100; 20; 100; 20; 100; 20; 100; 20];

% 2. Non-identical clocks
a0=0; % offset 100us 
b0=20; % skew 20ppm
a=500; % offset inittial value variation range: +/- 40 uniform distribution
b=30; % skew initial value variation raange: +/- 20 uniform distribution
iniOffset=a0+a*2*(rand(1,nNode)-0.5); stdOffset=sqrt((2*a)^2/12); % std for uniform dist.
iniSkew=b0+b*2*(rand(1,nNode)-0.5); stdSkew=sqrt((2*b)^2/12); %  sqrt((max-min)^2/12)std for uniform dist.
x0=10^(-6)*reshape([iniOffset;iniSkew],[],1);
% x0=10^(-6)*[10:10:2*nNode*10];
% illusrate the distributation of initial values 
% figure('name','distribution of inital values'); 
% subplot(2,1,1);hold on; plot(x0(1:2:end-1),'o-r'); title('Offset Initial Values'); grid on;ylabel('sec'); xlabel('node ID');
% subplot(2,1,2);hold on; plot(x0(2:2:end),'o-r'); title('Skew Initial Values');grid on;ylabel('ppm-like');xlabel('node ID');

% x0(2)=0; x0(1)=0; % set node 1 as ideal pefect clock
% x0=10^(-6)*[99.36; 19.25; 99.3; 19.3; 98.7; 18.5; 99.8; 19.8; 99.36; 19.5; ...
%             99.3; 19.3; 98.7; 18.5; 99.8; 19.8; 99.36;19.5; 99.3;  19.3;
%             98.7; 18.5; 99.8; 19.8];
fprintf("Clock initial values:\r");
fprintf("    Offset %d +/- %d us (std=%d us) \r", a0, a, stdOffset);
fprintf("    Skew  %d +/- %d ppm  (std=%d ppm)\r ", b0, b,stdSkew);

% internal variables for clock state, output and sync errors
% %noKalman
y=zeros(2*nNode,szsim); % output, a matrix for outputs at all sz simulation steps
xx=zeros(2*nNode,szsim); % state, a matrix for states at all sz simulation steps

yerr=zeros(2*nNode,szsim); % output errors

 indTheta=1:2:nNode*2-1;  % row indices for theta
 indSkew=2:2:nNode*2; % row indices for skew

% for root clock (node 1) failure
% disp("node 1 failure. Disconnected from network and ignore node 1 in statistis Ybar and std, etc");
% indTheta=3:2:nNode*2-1;  % row indices for theta
% indSkew=4:2:nNode*2; % row indices for skew


% xx(:,1)=x0_noisefree(1:2*nNode);
xx(:,1)=x0(1:2*nNode);
y(:,1)=xx(:,1)+measNoisev(:,1);
xx(:,1)=x0;
y(:,1)=x0;

% initial errors for interation from k=2  
yk=y(:,1);
Ybar(:,1)=[mean(yk(indTheta)');mean(yk(indSkew)')];
% yerr=[]; % clear yerr in case yerr has defined in previous simulation
yerr(:,1)=yk-kron(ones(nNode,1),Ybar(:,1));%误差向量（输出-输出均值）胡equ29  ??
   
%% Simulaiton Configuration 2d: Revise the node's role and topology
% update the nNode and nEdge, as they may be changed due to spannning tree
% nNode=numnodes(netG);% number of nodes
% nEdge=numedges(netG); % number of edges

% Revise the network to fit the topology needs
% Set the non-servo clock, reference clock and initial states
fprintf("Network size: %d nodes. Topology: L=", nNode);
% disp(L);
listnonServoClk=[];
listRefClk=[];
while true
    nID=input('Any clock node you want to change to non-sevo clock? \n     Type the node ID (0 for nothing to change): ');
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

% any changes to the network and clock configuraiton 
fprintf("Network created with : %d nodes.%d edges\n", nNode,nEdge);
%%  Simulating the Networked Synchronization controller 
disp('Now start simulation ');
hfigsim=figure('Name','Simulation Animation'); 
    cm = colormap('Lines'); 
    subplot(2,2,1);    title('offset \theta of all nodes & avg');
    subplot(2,2,2);    title('skew \gamma of all nodes & avg');
    subplot(2,2,3);    title('errors of \theta wrt.the average');
    % ax = axes; ax.ColorOrder = cm;
    subplot(2,2,4);    title('errors of \gamma wrt.the average');
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
    warning("     Unstable eigenvalue of the networked closed loop system \n");
else
    disp("     Good,stable system (two eigenvalue ==1, others <1)");
end
strK=num2str(K);
fprintf("     K=[%s ;\n        %s]\n", strK(1,:), strK(2,:));
fprintf("    simulaiton is in process");

% Please note, B1*L1 shoulde be the same as BK1
% sz=500;

for k = 2:szsim

%   % x-based simulation, State and output updates, Hu2019, eq.26, 27
    U=L1*y(:,k-1); % get the output differece with neighbours
    % U=L1*zeros(2*nNode,1); % free running clock, no regulation
    xx(:,k)=A1*xx(:,k-1)-B1*U+procNoisew(:,k-1); % state updates
    y(:,k)=xx(:,k)+measNoisev(:,k); % output updates

%     % y-base sys, noise-free      
%     y1(:,k)=A1*y1(:,k-1)-BK1*y1(:,k-1);%w输出Y值 胡 eq 28
%     xx1(:,k)=y1(:,k);%状态X值
    
    % collect the synchronization error for results analysis
    % Selet which reuslts to be used: 2 for x-based simulaiton, 1 for
    % y-based

    yk=y(:,k);
    Ybar(:,k)=[mean(yk(indTheta)');mean(yk(indSkew)')];
    Ystd(:,k)=[std(yk(indTheta)');std(yk(indSkew)')];
    yerr(:,k)=yk-kron(ones(nNode,1),Ybar(:,k));%误差向量（输出-输出均值）胡equ29  ??
   

   % for animaltion during simulation
    if (k<10) || (k<100 && (mod(k,10)==0)) ||(k>100 && mod(k,50)==0)
           figure(hfigsim);
           fprintf(".");

           % fprintf('k=%d\n',k);
        subplot(2,2,1);   drawTrajectory(xx(indTheta,1:k))
        hold on;  plot(1:k,Ybar(1,1:k), '-k', 'MarkerSize', 5, 'LineWidth', 2);
        subplot(2,2,2);   drawTrajectory(xx(indSkew,1:k))
        hold on;  plot(1:k,Ybar(2,1:k), '-sk', 'MarkerSize', 5,'LineWidth', 2);
        subplot(2,2,3);
        hold on; cla; plot(1:k,yerr(indTheta,1:k), '-o');
        subplot(2,2,4);
        hold on; cla; plot(1:k,yerr(indSkew,1:k), '-o','LineWidth', 1);
        title('errors of \gamma wrt.the average');
        xlabel('time (s)');
        drawnow
    end
    
end

fprintf('\n Simulation Ends\n');
%% plot results
simNotes=sprintf('Clk %d, ideal ref clk',1); xpos=1;nlayer=5;

if exist('Kname','var')    % for plot 
    Kname='K=xx'; NoiseSeq=0; 
end
strConfig=['(' Kname ', noise' num2str(NoiseSeq) ')'];
figure('name',['Output offset and offset error' strConfig]);
 subplot(2,2,1); plot(y(indTheta,:)','.-'); title('Offset output');
 hold on;  plot(1:k,Ybar(1,:), '-k', 'MarkerSize', 5, 'LineWidth', 2);
 subplot(2,2,2); plot(y(indSkew,:)','.-'); title('Skew output');
 hold on;  plot(1:k,Ybar(2,:), '-sk', 'MarkerSize', 1,'LineWidth', 2);
 subplot(2,2,3); plot(yerr(indTheta,1:100)', '.-'); title('Offset Errors');
 subplot(2,2,4); plot(yerr(indSkew,1:100)', '.-');title('Skew Errors');
  xlabel('time (s)');

  

figure('Name',['Net Graph & Offset/Skew variance ' strConfig]); 
k=length(Ystd);
subplot(1,2,1);
if issymmetric(L)
   plot(netG); %undirected graph
   % plot(netG,'LineWidth',2);
   % set(gca,'ydir','reverse'); % use ydir, xdir to adjust the figure's dir
   % set(gca,'xdir','reverse'); 
else
   % idrected graph. transpose the Laplacian for correct visulzation of
   % link's directiorn
   d=diag(L);     % D is the degree matrix 
    AdjM=diag(diag(L))-L; 
   plot(digraph(AdjM));
end
text(xpos-0.2,nlayer+0.2,simNotes);
% title(sprintf('Topology: %d nodes, %d edges (ref clk = [%s])',nNode, nEdge, num2str(listnonServoClk)));
if (length(listRefClk)==0)
   title(sprintf('Topology: %d nodes, %d edges (no refence clk)',nNode, nEdge)); 
else
   title( "Topology: "+nNode+" nodes  "+nEdge+ "edges  "+ num2str(length(listRefClk)+ " Refence clock = [" + num2str(listRefClk))+"]");
end 

%offset的方差以及滑动平均下的方差
subplot(1,2,2);hold on;
plot(1:k,Ystd, 'MarkerSize', 5, 'LineWidth', 1);
hold on;
title(['std of offset and skew (' Kname ', n' num2str(NoiseSeq) ')']);
mov50Offset=movmean(Ystd(1,:),50); 
mov50Skew=movmean(Ystd(2,:),50); 
plot(1:k,mov50Offset,'b','LineWidth', 2);%线性搜索的 
plot(1:k,mov50Skew,'r','LineWidth', 2)%线性搜索的

xlabel('time (s)');
legend('std of offset','std of skew','moving avg of offset std','moving avg of skew std');
grid on;
% axes('position',[.10 .5 .2 .25]);%第二个图摆放的横向位置、纵向位置、宽度、高度
% box on
% M1=movmean(Ystd(1,:),10);plot(1:k,M1,'g')%传统的
% M1=movmean(Ystd(1,:),10);plot(1:k,M1,'b')%线性搜索的
% axes('position',[.7 .5 .2 .25]);%第二个图摆放的横向位置、纵向位置、宽度、高度
% box on
% M2=movmean(Ystd(2,:),10);plot(1:k,M2,'r')%传统的
% M2=movmean(Ystd(2,:),10);plot(1:k,M2,'y')%线性搜索的



  
fprintf("Simulation Done!  %d nodes, %d sync steps\n",nNode,szsim);
 disp("     " + num2str(length(listnonServoClk)+ " non-servo clock = [" + num2str(listnonServoClk))+"]");
 disp("     " + num2str(length(listRefClk)+ " Refence clock = [" + num2str(listRefClk))+"]");
% fprintf("   mean value of all nodes at end: [offsetba=%d skewbar=%d]\n",Ybar(1,szsim), Ybar(2,szsim));
% fprintf("   with standard deviation [%d, %d]\n", Ystd(1,k), Ystd(2,k));
fprintf("   maximum std in t=[300,800]: [offset=%d skew=%d]\n",max(Ystd(1,[300:szsim])), max(Ystd(2,[300:szsim])));
fprintf("   average std in t=[300,800]: [offset=%d, skew=%d]\n", mean(Ystd(1,[300:szsim])), mean(Ystd(2,[300:szsim])));


xvar=var(xx(:,[300:szsim])');
xthetavar=[1:nNode; xvar(1:2:2*nNode-1)];
xskewvar=[1:nNode; xvar(2:2:2*nNode)];

% xvarsorted=sortrows(xthetavar',[2]); % sort in ascending order based on the elements in the 2 column.
% disp("node (in ascending theta var) ="); disp(xvarsorted(:,1)');
% disp("                    theta var ="); disp(xvarsorted(:,2)');

disp("Simulaiton Results Statistics End");
return;
%%  plot results for visualization
figure('name','Sychronization offset error');
hold on; cla; plot(1:k,yerr(indTheta,1:k), '-.');

      
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
      

        
        


        
   
