figure('name','Std of offsets');
load Result_Kold_huan20.mat
k=length(Ystd);
hold on;
plot(1:k,Ystd(1,:)*1e6,'b', 'MarkerSize', 5, 'LineWidth', 1);
plot(1:k,mov50Offset*1e6,'-.b','LineWidth', 2)
hmain = gca;

xlabel('Time k (sec)');
ylabel('std of { \theta_i (k) } (\mus)');
title('Standard Deviation'); 
legend('conventional method K0','50s moving average for K0','linear search optimization K1','50s moveing aeverage for K1');
grid on;

hax=axes('position',[.38 .3 .5 .25]);
zoominRange=[200:k];
plot(zoominRange,Ystd(1,zoominRange)*1e6,'b', 'MarkerSize', 5, 'LineWidth', 1);
hold on;
plot(zoominRange,mov50Offset(zoominRange)*1e6,'-.b','LineWidth', 2)
box on;
title("Zoomed in [300,800] sec");
ylabel("\mus");
xlabel("Time (sec");


%% 
figure('name','Std of skew');
load Result_Kold_huan20.mat
k=length(Ystd);
hold on;
plot(1:k,Ystd(2,:)*1e6,'r', 'MarkerSize', 5, 'LineWidth', 1);
plot(1:k,mov50Skew*1e6,'-.r','LineWidth', 2)
hmain = gca;

xlabel('Time k (sec)');
% ylabel('std of { \theta_i (k) } (\mus)');
ylabel('std of { \gamma_i (k) } (ppm)');
title('Standard Deviation'); 
legend('conventional method K0','50s moving average for K0');
grid on;

hax=axes('position',[.38 .3 .5 .25]);
zoominRange=[200:k];
plot(zoominRange,Ystd(2,zoominRange)*1e6,'r', 'MarkerSize', 5, 'LineWidth', 1);
hold on;
plot(zoominRange,mov50Offset(zoominRange)*1e6,'-.r','LineWidth', 2)
box on;


xlabel('Time k (sec)');
% ylabel('std of { \theta_i (k) } (\mus)');
ylabel('std of { \gamma_i (k) } (ppm)');
title('Standard Deviation'); 
legend('conventional method K0','50s moving average for K0');
grid on;






%% plot results of linear search 
     figure('name','K');
     plot(allK(:,1,1),'.-b'); hold on;
     plot(allK(:,1,2),'.-r');
     plot(allK(:,2,1),'o-b');
     plot(allK(:,2,2),'o-r');
     title('gain matrix K (blue: k11,k21, red k12,k22)');
     xlabel('training episode');
     
     figure('name','Ycopt');
     plot(J1,'-ob','LIneWidth',2); title('J1=Ycopt'); xlabel('main iteration i');
%       plot(1:2:14,'-ob','LIneWidth',2); title('J1=Ycopt'); xlabel('main iteration i');
     % pause(0.2);
     
 