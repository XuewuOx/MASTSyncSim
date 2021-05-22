% % ***************************************************************************
% * Created on:  31/03/2021
% * Modified on: 13/05/2021
% * 
% * Author:      Yan Zong, Xuewu Dai
% *
% * Descriptoin: the code in this branch is used to simulate the network
% *              under the static controller, the corresponding paper is
% *              submitted to IEEE Internet of Things J. Special Issue for 
% *              peer review.
% *
% * IDE:         MATLAB R2020b
% % ***********************************************************************
clear all;
close all;
clc;
%% Simulaiton Configuration 1: Network Topology
disp("Clock Synchronization Simulation");

szsim=800; % simulation time

load('OriginalDirectedGraphwith50Nodes&60Edges.mat','L');
load('OriginalDirectedGraphwith50Nodes&60Edges.mat','nNode');

fprintf("The number of tree spinning network is %d \r", nNode);

[NetTree]=genTreeNet(L);

% directed graph
d=diag(NetTree);  % d is the dialog vector of matrix L 
A=diag(d)-NetTree;
netG=digraph(A,'omitselfloops');
figure('name', 'Network Topology'); plot(netG,'LineWidth',2);
title("a directed  graph ");

% disp(NetTree);
%% Simulaiton Configuration 2a: Clock and Packet-exchange Delay Noises (Q)
% According to the 3Sigma rule of Gaussian distribution 99.73% of the noise 
% will be in the range of [mu-3*sigma, mu+3*sigma] This suggests nearly all 
% values are taken to lie within three standard deviations of the mean, and 
% thus it is empirically useful to treat 99.7% probability as near certainty.
% Therefore, the noise range is -/+ 3*sigma
%
% A rule of thumb, according to which, in certain problems in probability 
% theory and mathematical statistics, an event is considered to be 
% practically impossible if it lies in the region of values of the normal 
% distribution of a random variable at a distance from its mathematical 
% expectation of more than three times the standard deviation.
% 
% In the mathematical version, 3 signms rule can be expressed as the
% following: 
% If abs(obs-mean(obs))>3*sd(obs) ===> you've got something unusual.
%

sigma1sqr=10^-12; % variance of offset noise, from Clock A of Giorgi2011
sigma2sqr=10^-16; % variance of skew noise, from Clock A of Giorgi2011
sigma3sqr=(4*10^(-6))^2; % variance of measurement noise, from Zong2019c

mu=[0 0];
R=[sigma1sqr 0;0 sigma2sqr]; % covariance of process noise (NO USE)
Q=[sigma3sqr sigma3sqr;sigma3sqr 2*sigma3sqr]; % covariance of measurement noise

fprintf("noise levels:\r");
fprintf("    offset process noise std = %d\r", sqrt(sigma1sqr));
fprintf("    skew process noise std = %d\r", sqrt(sigma2sqr));
fprintf("    offset measurement noise std = %d\r", sqrt(sigma3sqr));

% for the k-th node, the noises is a 2-by-sz matrix stated in  
% odd row of procNoise are (1) process noses of theta  (1-by-sz) 
% even row of procNoise are (2) process noise of gamma (1-by-sz)
% odd row of measNoise are (1) measurement noise of theta (1-by-sz)
% odd row of measNoise are (2) measurement noise of gamma (1-by-sz)

procNoise=[];
measNoise=[];

for k=1:nNode
    procNoise((k-1)*2+1:(k-1)*2+2,:)=[sqrt(sigma1sqr)*randn(1,szsim);sqrt(sigma2sqr)*randn(1,szsim)];
    measNoise((k-1)*2+1:(k-1)*2+2,:)=mvnrnd(mu,Q,szsim)'; 
end
%% Simulaiton Configuration 2b: Controller Design  
% static control gain is obtained by using the LMI technique
% the control gain from LMI of Chang2014
K=[-0.0021 0; 
    0.0001 -0.0026];

% the control gain from IEEE Trans. Cybern.
K=[0.5 0;
    0 0.025];

format long   
fprintf("Static controller gain K is given by using LMI:\n"); disp(K);
format short
%% Simulaiton Configuration 2c: Clock & Networked State Initialisation
% initialising the clock offset and skew (between 0 and 50ppm)
a0=600000; % initial offset is 600ms 
b0=25; % initial skew is 25ppm
a=200000; % initial offset variation range: +/- 200 ms uniform distribution
b=25;   % initial skew variation raange: +/- 25 uniform distribution
theta0=a0+a*2*(rand(1,nNode)-0.5); 
stdOffset=sqrt((2*a)^2/12); % std. dev. for uniform dist. sqrt((max-min)^2/12)
gamma0=b0+b*2*(rand(1,nNode)-0.5); 
stdSkew=sqrt((2*b)^2/12); % std. dev. for uniform dist.
x0=10^(-6)*reshape([theta0;gamma0],[],1); % reshape matrix to have specified 
                                              % number (i.e. 1) of columns,
                                              % the unit is "us"
fprintf("initial system state x[0]:\r");
fprintf("    clock offset %d +/- %d us, & std=%d us \r", a0, a, stdOffset);
fprintf("    clock skew %d +/- %d ppm, & std=%d ppm \r", b0, b, stdSkew);

T=1; % clock synchronization interval
A=[1 T;0 1];
B=[1 0;0 1];

% internal variables for clock state, output and sync errors
x=zeros(2*nNode,szsim); % state, a matrix for states at all sz simulation steps
y=zeros(2*nNode,szsim); % output, a matrix for outputs at all sz simulation steps
yerr=zeros(2*nNode,szsim); % synchronisation error

indTheta=1:2:nNode*2-1;  % row index for theta
indSkew=2:2:nNode*2; % row index for skew

% for first row
x(:,1)=x0(1:2*nNode); % system state
y(:,1)=x(:,1)+measNoise(:,1); % system output
yy(:,1)=x(:,1); % for performance evaluation

% initial errors for interation from k=2  
yk=yy(:,1);
Ybar(:,1)=[mean(yk(indTheta)');mean(yk(indSkew)')];
yerr(:,1)=yy; % synchronisation error (actually is the clock offset)
   
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
    refClk=input(sprintf('Set node %d as a reference clock (theta=0, gamma=0, procNoise=0) (Y) or just non-servo clock (N)? ',...
            nID),'s');
    if refClk=='y' || refClk=='Y'
        x0((nID-1)*2+1:(nID-1)*2+2)=0; % initial offset and gamme to zero
        xx2(:,1)=x0;
        y2(:,1)=x0;
        procNoise([(nID-1)*2+1:(nID-1)*2+2],:)=0; % clear process noises w
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
    subplot(2,2,4);    title('errors of \gamma wrt.the average');

A1=kron(eye(nNode),A); % Kronecker product(克罗内克积)
BK=B*K;
BK1=kron(L,BK);    
B1=kron(eye(nNode),B*K);
L1=kron(L,eye(2)); % B1*L1 is the same as BK1

if chkEigAc(A,B,K,L)==false
    % the close form of the NCS has unstable eigenvalue
    warning("     Unstable eigenvalue of the networked closed loop system \n");
else
    disp("     Good,stable system (two eigenvalue ==1, others <1)");
end
strK=num2str(K);
fprintf("     K=[%s ;\n        %s]\n", strK(1,:), strK(2,:));
fprintf("    simulaiton is in process");

for k = 2:szsim

    % state and output updates, see eq.26, 27 in Hu2019
    U=L1*y(:,k-1); % get the output differece with neighbours
    x(:,k)=A1*x(:,k-1)-B1*U+procNoise(:,k-1); % state updates
    y(:,k)=x(:,k)+measNoise(:,k); % output updates

    % collect the synchronisation error for result analysis
    yy=x;
    yk=y(:,k);
    Ybar(:,k)=[mean(yk(indTheta)');mean(yk(indSkew)')];
    Ystd(:,k)=[std(yk(indTheta)');std(yk(indSkew)')];
    yerr(:,k)=yy(:,k);
   
   % for animaltion during simulation
    if (k<10) || (k<100 && (mod(k,10)==0)) ||(k>100 && mod(k,50)==0)
           figure(hfigsim);
           fprintf(".");

           % fprintf('k=%d\n',k);
        subplot(2,2,1);   drawTrajectory(x(indTheta,1:k))
        hold on;  plot(1:k,Ybar(1,1:k), '-k', 'MarkerSize', 5, 'LineWidth', 2);
        subplot(2,2,2);   drawTrajectory(x(indSkew,1:k))
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


xvar=var(x(:,[300:szsim])');
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
            plot([1:szsim],x(kk,:));hold on
        end
      %    legend('节点1 offset状态值','节点2 offset状态值','节点3 offset状态值','节点4 offset状态值','节点5 offset状态值','节点6 offset状态值','节点7 offset状态值','节点8 offset状态值','节点9 offset状态值','节点10 offset状态值','节点11 offset状态值','节点12 offset状态值');
        xlabel('n');ylabel('状态值');
        
        subplot(2,2,2)
        for kk=2:2:2*nNode
            plot([1:szsim],x(kk,:));hold on
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
            var_xx(kk)=var(x(kk,:));
        end
      mean_xx=mean(x,2);
      

        
        


        
   
