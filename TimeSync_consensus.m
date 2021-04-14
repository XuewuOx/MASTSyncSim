% % ***********************************************************************
% * Descriptoin: the code in this branch is used to simulate the network
%               under the PI controller, the corresponding paper is
%               submitted to IEEE Trans. Ind. Informat. for peer review
% * Limitation: the code can run on the MATLAB R2020b
% % ***********************************************************************
clear all;
close all;
clc;
%% Simulaiton Configuration 1: Network Topology
disp("Clock Synchronization Simulation");
sz=600; % simulation time
% linear network topology specified by the Laplacian matrix L
% bi-directed coupling network
% L=[1 -1 0 0 0 0 0   ;
%   -1 2 -1 0 0 0 0   ;
%   0 -1 2 -1 0 0 0   ;
%   0 0 -1 2 -1 0 0   ;
%   0 0 0 -1 2 -1 0   ;
%   0 0 0 0 -1 2 -1   ;
%   0 0 0 0 0 -1 1   ;];
% directed coupling network
L=[1 -1  0  0  0  0  0  0  0;
  -1  1  0  0  0  0  0  0  0;
   0 -1  1  0  0  0  0  0  0;
   0  0 -1  1  0  0  0  0  0;
   0  0  0 -1  1  0  0  0  0;
   0  0  0  0 -1  1  0  0  0;
   0  0  0  0  0 -1  1  0  0;
   0  0  0  0  0  0 -1  1  0;
   0  0  0  0  0  0  0 -1  1;];
nNode=length(L);    % number of nodes in a network
L(1,:)=zeros(1,nNode); % set the node[1] as the reference node, note that 
                       % since there is no index 0 in the matlab, node[1] 
                       % is seen as node[0]
disp(sprintf("Network size: %d nodes. Topology: L=", nNode));
disp(L);

%% Simulaiton Configuration 2: Clock and controller
% omega_theta is non-zero item, but the omega_gamma is zero. Since the 
% the item of omega_gamma corresponds to integeral controller in my paper
A=[1 1;0 1];
B=[1 0;0 1];
alpha=0.5;
beta=0.025;
K=[alpha 0;beta 0]; % the control gain is the same as that in the experiments 
fprintf("Control Gain K="); disp(K);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) initialising the clock offset and skew (less than 10ppm)
a0=600000; % initial offset is 100us 
b0=5;   % initial skew is 5ppm
a=200000;   % initial offset variation range: +/- 100 uniform distribution
b=5;   % initial skew variation raange: +/- 5 uniform distribution
theta0=a0+a*2*(rand(1,nNode)-0.5);  % initial clock offset
gamma0=b0+b*2*(rand(1,nNode)-0.5);  % initial clock skew
w0=0+0*2*(rand(1,nNode)-0.5);   % integral controller
% reshape matrix to have specified number of columns
x0=10^(-6)*reshape([theta0;w0],[],1); 
x0=x0';
fprintf("initial system state x[0]:\r");
fprintf("    clock offset %d +/- %d us \r", a0, a);
fprintf("    integral controller %d +/- %d \r ", 0, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) generating the process and measurement noise
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

% processing noise corresponds to the o[k] of equ(21) in the paper,
% measurement noise is the standard deviation of packet exchaneg delay in
% the paper, since we use the feedforward control to compensate for the
% effects of the packet exchange delay, thus we can use the standard
% deviation of the packet exchange delay to denote the measurement noise.

% for the k-th node, the noises is a 2-by-sz matrix stated in  
% odd row of procNoise are (1) process noses of theta  (1-by-sz) 过程噪声中thet的随机值
% even row of procNoise are (2) process noise of gamma (1-by-sz)  过程噪声中integralcontroller的随机值
% odd row of measNoise are (1) measurement noise of theta (1-by-sz) 测量噪声中theta 的随机值(生成多维正态数据）
% odd row of measNoise are (2) measurement noise of gamma (1-by-sz)

% yan: from eq (21), we can see that there is a value for compensating for
% the joint effects due to the desynchronisation and clock skew, we cannot 
% model this procedure in this code, since the scheduling issue in the
% network is modelled as the reference input convergence problem in the
% control theory. For the processing noise, this simulation only consider
% processing delay eta, and offset noise (i.e. the effects of skew and 
% phaise noise)

% the following two noise sigmas are for the processing noise
sigma1sqr=10^-12 + 4*10^-12;   % variance of offset noise, see (21), 
                               % adoped from Giorgi2011, this value is the
                               % sum of phase noise and processing noise
                               % eta
sigma2sqr=0;    % no process noise on the integral controller 

% 0.3^-6 is the standard deviation of the packet exchange delay (i.e.
% measurement noise. 
sigma3sqr=(4*10^(-6))^2;  % 1.6*10^-11
mu=[0 0];
sigma=sigma3sqr*[1 0;0 1];

fprintf("noise levels:\r");
fprintf("    offset process noise std = %d\r", sqrt(sigma1sqr));
fprintf("    integral controller process noise std = %d\r", sqrt(sigma2sqr));
fprintf("    offset measurement noise std = %d\r", sqrt(sigma3sqr));
                            
procNoise=[];
measNoise=[];
for k=1:nNode
procNoise((k-1)*2+1:(k-1)*2+2,:)=[sqrt(sigma1sqr)*randn(1,sz)+(gamma0(k)/1000000+311/1000000)*ones(1,sz);sqrt(sigma2sqr)*randn(1,sz)];
measNoise((k-1)*2+1:(k-1)*2+2,:)=mvnrnd(mu,sigma,sz)'; 
measNoise((k-1)*2+2,:)=measNoise((k-1)*2+1,:);  % the measurement noise on 
                                                % theta and integral controller 
                                                % should be same
end

% internal variables for clock state, output and sync error
y=zeros(2*nNode,sz);   % output, a matrix for output at all sz simulation steps
x=zeros(2*nNode,sz);  % state, a matrix for state at all sz simulation steps

indTheta=1:2:nNode*2-1; % row indices for theta (e.g. 1,3,5,...)
indw=2:2:nNode*2;    % row indices for integral controller (e.g. 2,4,6,...)
    
% for first row
x(:,1)=x0(1:2*nNode);   % system state
y(:,1)=x(:,1)+measNoise(:,1);  % system output
yy(:,1)=x(:,1); % for performance evaluation

x0(1:2)=0; % initial offset and gamme to zero (first column, first two rows)
x(:,1)=x0; % first column of x is assigned as first row of x0
y(:,1)=x0; % first column of y is assigned as first row of x0
procNoise([1:2],:)=0; % clear process noises (clear first two rows)

% initial errors for interation from k=2  
yk=yy(:,1);
Ybar(:,1)=[mean(yk(indTheta)');mean(yk(indw)')];
yerr(:,1)=yy;   % synchronisation error (actually is the clock offset)
%% Simulating the Networked Synchronization controller 
figure('Name','Simulation Animation'); 
    cm = colormap('Lines'); 
    subplot(2,2,1);
    title('offset \theta of all nodes & avg');
    subplot(2,2,2);
    title('integral controller w of all nodes');
    subplot(2,2,3);
    title('time synchronisation error (i.e. \theta)');
%     subplot(2,2,4);
%     title('errors of \gamma wrt.the average ToDo');

A1=kron(eye(nNode),A);  % Kronecker product(克罗内克积)
BK=B*K;  % variable D for the feedback gain matrix K in Hu2019
BK1=kron(L,BK);  
B1=kron(eye(nNode),B*K);
L1=kron(L,eye(2));  % B1*L1 is the same as BK1

for k = 2:sz

    % state and output updates, see eq.26, 27 in Hu2019
    U=L1*y(:,k-1); % get the output differece with neighbours
    x(:,k)=A1*x(:,k-1)-B1*U+procNoise(:,k-1); % state updates
    y(:,k)=x(:,k)+measNoise(:,k); % output updates
    
    % collect the synchronization error for results analysis
    yy=x;
    yk=yy(:,k);
    Ybar(:,k)=[mean(yk(indTheta)');mean(yk(indw)')];
    Ystd(:,k)=[std(yk(indTheta)');std(yk(indw)')];    
    yerr(:,k)=yy(:,k);
   
    subplot(2,2,1);
    drawTrajectory(x(indTheta,1:k))
    hold on;  plot(1:k,Ybar(1,1:k), '-k', 'MarkerSize', 5, 'LineWidth', 2);

    subplot(2,2,2);
    drawTrajectory(x(indw,1:k))
    hold on;  plot(1:k,Ybar(2,1:k), '-sk', 'MarkerSize', 5,'LineWidth', 2);
    
    subplot(2,2,3);
    hold on; cla; plot(1:k,yerr(indTheta,1:k), '-o');
    
%     subplot(2,2,4);
%     hold on; cla; plot(1:k,yerr(indw,1:k), '-o','LineWidth', 1);
   
    if (k<10) || (k<100 && (mod(k,10)==0)) ||(k>100 && mod(k,50)==0)
       drawnow
    end

end
fprintf("Simulation Done!  %d nodes, %d sync steps\n",nNode,sz);
fprintf("    mean value of all nodes at end: [offset_ba=%d]\n",Ybar(1,sz));
fprintf("    with standard deviation [%d]\n", Ystd(1,k));

cd 'C:\Users\yan\Documents\Projects\ctrl_anlys\LinearSearchIteration';
save('simulation_all_parameters');
save('ts_precision','yerr');