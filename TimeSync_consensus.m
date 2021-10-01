% % ***************************************************************************
% * Created on:  31/03/2021
% * Modified on: 13/05/2021
% * 
% * Author:      Yan Zong, Xuewu Dai
% *
% * Descriptoin: the code in this branch is used to compare TPSN, dynamic 
% *              controller, moving average and PISync, the corresponding 
% *              paper is submitted to IEEE TII for peer review.              
% *
% * IDE:         MATLAB R2020b
% % ***********************************************************************
clear all;
close all;
clc;

% selecting the different clock offset and skew processing methods:
% 1 -- enable; 0 -- disable
PISYNC = 1; 
DYNCTRL = 0;
MOVAVG = 0;
TPSN = 0;
%% Simulaiton Configuration 1: Network Topology
disp("Clock Synchronisation Simulation");

szsim=1000; % simulation time

load('OriginalDirectedGraphwith50Nodes&60Edges.mat','L');
load('OriginalDirectedGraphwith50Nodes&60Edges.mat','nNode');

fprintf("The number of tree spinning network is %d \r", nNode);

[NetTree]=genTreeNet(L);

% directed graph
d=diag(NetTree);  % d is the dialog vector of matrix L 
AA=diag(d)-NetTree;
netG=digraph(AA,'omitselfloops');
% figure('name', 'Network Topology'); 
GG = plot(netG,'Layout','layered','Direction','up','LineStyle','--','NodeFontName','Times New Roman','NodeFontSize',9.5,'Interpreter','latex'); 
% by default 'Linewidth' is 0.5, and 'MarkerSize' is 4 

highlight(GG,[1],'NodeColor',[0.6350 0.0780 0.1840], 'MarkerSize', 5.5) 
NodeIndex=1:1:50;
for ii=1:50
    NewNodeIndexTemp = num2str(NodeIndex(ii)-1);
    NewNodeIndex(ii) = cellstr(NewNodeIndexTemp);
end 
labelnode(GG,NodeIndex,NewNodeIndex)

set(gca, 'FontSize', 11, 'FontName', 'Times New Roman')  % axis configuration

% set(gca, 'xlim',[0 0.7], 'XTick',0:0.1:0.7, 'LineWidth', 1.5); % axis configuration
% set(gca,'ylim',[-40 80], 'YTick',-40:20:80); % axis configuration

% xlabel('Xx', 'Interpreter','latex', 'FontSize', 11, 'FontName', 'Times New Roman');
% ylabel('Yy', 'Interpreter','latex', 'FontSize', 11, 'FontName', 'Times New Roman');
% title("Network Topology");

set(gcf,'color','w');
set(gcf,'renderer','Painters');
% print -depsc -tiff -r600 -painters fig2.eps;
% saveas(gcf,'fig2.tif');

disp(NetTree);
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
sigma2sqr=10^-12; % variance of skew noise, from Clock A of Giorgi2011
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the control gain from Yildirim2018, PISync
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (PISYNC == 1 | TPSN == 1)
    alpha = 1;
    
    if PISYNC == 1
        beta = 0; % see equ (9) of Yildirim2018
        initial_beta = 1/32768000;
    elseif TPSN == 1
        beta = 1;
    end 
    
    K = [alpha 0; 0 beta];

    format long   
    fprintf("Static controller gain K is given by using LMI:\n"); disp(K);
    format short
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamic controller gain obtained by using the LMI technique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (DYNCTRL == 1 | MOVAVG == 1)
    % K = [A_K B_K; C_K D_K]
    if DYNCTRL == 1
        K=-[-0.039511041468828 0 -0.000000000025197 0;
            0 -0.039511041572130 0 0.000000000010100;
           -0.020919912004335 0 -0.775989213038119 0;
           0 -0.004456899766823 0 -0.761027326480981]; 
        Gamma=5.3330;
    elseif MOVAVG == 1
        K=[0   0   0   0; % The configurarions of moving average is from Tian2021
           0   0.5 0   0.5;
           0   0   0.5   0;
           0   0.5 0   0.5]                
    end
    
    A_K = K(1:2, 1:2);
    B_K = K(1:2, 3:4);
    C_K = K(3:4, 1:2);
    D_K = K(3:4, 3:4);    
        
    format long   
    fprintf("Dynamic controller gain A_K, B_K, C_K, D_K are given by using LMI:\n"); 
    fprintf("A_K = \n"); disp(A_K);
    fprintf("B_K = \n"); disp(B_K);
    fprintf("C_K = \n"); disp(C_K);
    fprintf("D_K = \n"); disp(D_K);
    
    if DYNCTRL == 1
        fprintf("The H_infty performance gamma = \n"); disp(Gamma);
    end
    
    format short
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
Beta=zeros(nNode,szsim); % history of tunned betas

indTheta=1:2:nNode*2-1;  % row index for theta
indSkew=2:2:nNode*2; % row index for skew

% for first row
x(:,1)=x0(1:2*nNode); % system state
y(:,1)=x(:,1)+measNoise(:,1); % system output
yy(:,1)=x(:,1); % for performance evaluation

% setting node 1 as a reference node
x0(1:2)=0; % initial offset and gamme to zero (first column, first two rows)
x(:,1)=x0; % first column of x is assigned as first row of x0
y([1:2],:)=0; % first two rows of y is assigned as 0 (clear first two rows) 
yy(:,1)=x0; % first column of yy is assigned as first row of x0
procNoise([1:2],:)=0; % clear process noises (clear first two rows)
measNoise([1:2],:)=0; % clear process noises (clear first two rows)

% initial errors for interation from k=2  
yk=yy(:,1);
Ybar(:,1)=[mean(yk(indTheta)');mean(yk(indSkew)')];
yerr(:,1)=yy; % synchronisation error (actually is the clock offset)
%%  Simulating the Networked Synchronization controller 
disp('Now start simulation ');
hfigsim=figure('Name','Simulation Animation'); 
    cm = colormap('Lines'); 
    subplot(2,2,1);    title('offset \theta of all nodes & avg');
    subplot(2,2,2);    title('skew \gamma of all nodes & avg');
    subplot(2,2,3);    title('errors of \theta wrt.the average');
    subplot(2,2,4);    title('errors of \gamma wrt.the average');

A1=kron(eye(nNode),A); % Kronecker product(克罗内克积)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the control gain from Yildirim2018, PISync
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (PISYNC == 1 | TPSN == 1)
    % B1=kron(eye(nNode),B*K); % all the initial control gains are same
    alpha_list = alpha * ones([nNode, 1]);
    beta_list = beta * ones([nNode, 1]);
    B1 = zeros([2*nNode,2*nNode]);
    for i=1:nNode
        B1((i-1)*2+1, (i-1)*2+1) = alpha_list(i, 1);
        B1((i-1)*2+2, (i-1)*2+2) = beta_list(i, 1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamic controller gain obtained by using the LMI technique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (DYNCTRL == 1 | MOVAVG == 1)   
    x_K=zeros(2*nNode,szsim);
    A_K1=kron(eye(nNode),A_K);
    B_K1=kron(eye(nNode),B_K);
    C_K1=kron(eye(nNode),B*C_K);
    D_K1=kron(eye(nNode),B*D_K);
end 

NetTreeTemp=kron(NetTree,eye(2)); % B1*L1 is the same as BK1

% if chkEigAc(A,B,K,NetTree)==false
%     % the close form of the NCS has unstable eigenvalue
%     warning("     Unstable eigenvalue of the networked closed loop system \n");
% else
%     disp("     Good,stable system (two eigenvalue ==1, others <1)");
% end
strK=num2str(K);
fprintf("     K=[%s ;\n        %s]\n", strK(1,:), strK(2,:));
fprintf("    simulaiton is in process");

for k = 2:szsim

    % state and output updates, see eq.26, 27 in Hu2019
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    % dynamic controller gain obtained by using the LMI technique
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if (DYNCTRL == 1 | MOVAVG == 1)
        UTmp=NetTreeTemp*y(:,k-1); % get the output differece with neighbours
        x_K(:,k)=A_K1*x_K(:,k-1)+B_K1*UTmp; % x_F[k+1] = A_F * x_F[k] + B_F * y[k]    
        U = C_K1*x_K(:,k)+D_K1*UTmp; % u[k] = C_F * x_F[k] + D_F * y[k]
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % the control gain from Yildirim2018, PISync
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    if (PISYNC == 1 | TPSN == 1)
        UTmp=NetTreeTemp*y(:,k-1); % get the output differece with neighbours
        U = B1*UTmp;
    end
    
    % state and output updates, see eq.26, 27 in Hu2019      
    x(:,k)=A1*x(:,k-1)-U+procNoise(:,k-1); % state updates
    y(:,k)=x(:,k)+measNoise(:,k); % output updates
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % update the controlling gain by using the PISync protocol
    % alpha is always equal to one in PISync
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if PISYNC == 1
        if k>3
            y_product = (y(:,(k)) - y(:,k-1)) .* (y(:,(k-1)) - y(:,k-2));
        else                    
            if k == 2
                y_product = (y(:,(k)) - y(:,k-1)) .* (y(:,(k-1)) - zeros(2*nNode,1));    
            elseif k == 1
                y_product = zeros(2*nNode,1);
            end 
        end

        for i=2:nNode       

            if y((i-1)*2+1, k) > (b*2*2 / 32768)        
                beta_list(i, 1) = 0;
            elseif beta_list(i, 1) == 0
                beta_list(i, 1) = initial_beta;
            elseif y_product((i-1)*2+1, 1) > 0
                beta_list(i, 1) = beta_list(i, 1) * 2; % lambda^+ = 2
                beta_list(i, 1) = max(beta_list(i, 1), initial_beta);
            else
                beta_list(i, 1) = beta_list(i, 1) / 3; % lambda^- = 3            
            end 

            B1((i-1)*2+1, (i-1)*2+1) = alpha_list(i, 1);
            B1((i-1)*2+2, (i-1)*2+2) = beta_list(i, 1);
        end   
        Beta(:,k)=beta_list;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
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

%%
% plotting the offset and skew outputs
figure

offset_us=zeros(szsim, nNode);
for i=1:1:nNode
    offset_us(:, i) = 1000000 * yerr((i-1)*2+1,:)'; % change the unit to us
end 

subplot(2,2,1); 
for i=1:1:nNode
    plot(offset_us(:, i));
    hold on;
end 
set(gca, 'xlim',[0 60], 'XTick',0:10:60); % axis configuration
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel('$\theta_i[k]$ (s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
title('Clock offset at convergence state', 'Interpreter','latex');

subplot(2,2,2); 
for i=1:1:nNode
    plot(offset_us(:, i));
    hold on;
end 
set(gca, 'xlim',[szsim-600 szsim], 'XTick',szsim-600:100:szsim); % axis configuration
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel('$\theta_i[k]$ ($\mu$s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
title('Clock offset at steady state', 'Interpreter','latex');

skew_ppm=zeros(szsim, nNode);
for i=1:1:nNode
    skew_ppm(:, i) = 1000000 * yerr((i-1)*2+2,:)'; % change the unit to ppm 
end

subplot(2,2,3); 
for i=1:1:nNode
    plot(skew_ppm(:, i));
    hold on;
end 
set(gca, 'xlim',[0 60], 'XTick',0:10:60); % axis configuration
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel('$\gamma_i[k]$ (ppm)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
title('Clock Skew at convergence state', 'Interpreter','latex');

subplot(2,2,4); 
for i=1:1:nNode
    plot(skew_ppm(:, i));
    hold on;
end 
set(gca, 'xlim',[szsim-600 szsim], 'XTick',szsim-600:100:szsim); % axis configuration
xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel('$\gamma_i[k]$ (ppm)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
title('Clock Skew at steady state', 'Interpreter','latex');

%%
% plotting the tunning beta 
if PISYNC == 1
    figure
    for i=1:1:nNode
        plot(Beta(i,:));
        hold on;
    end
    xlabel('Time (s)', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
    ylabel('$\beta_i[k]$', 'Interpreter','latex', 'FontSize', 13, 'FontName', 'Times New Roman');
    title('Evolution of $\beta_i[k]$ in the simulations', 'Interpreter','latex');
end

%%
% saving the simulation results
if PISYNC == 1
    save('simulation_all_parameters_pisync');
    save('ts_precision_pisync','yerr');     
elseif DYNCTRL == 1
    save('simulation_all_parameters_dynctrl');
    save('ts_precision_dynctrl','yerr');     
elseif MOVAVG == 1
    save('simulation_all_parameters_movavg');
    save('ts_precision_movavg','yerr');        
elseif TPSN == 1
    save('simulation_all_parameters_tpsn');
    save('ts_precision_tpsn','yerr');    
end
