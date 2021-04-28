function [stableFlag, eigAnetc,numEigOne] = chkEigAc(A,B,K,L,debugplot)
%chkEigAc to check if the closed-loop network system is stable
%   by caculating the eigenvalue of the networked closed-loop system
%   state transition matrix  A-BK

%
if nargin==5
    DEBUGPLOT=(debugplot==1);
else
    DEBUGPLOT = false;
end
stableFlag=true; % default value
nA=size(A,1);
nNode=size(L,1);
Anet=kron(eye(nNode),A);
BK=B*K;  % feedback gain matrix K
BKnet=kron(L,BK);
eigAnetc=eig(Anet-BKnet);
eigR=abs(eigAnetc);
indOne=abs(eigR-1)<1e-8;
numEigOne=sum(indOne);
fprintf("\n  tor the system (A,B,K,L): number (multiplicity) of one eigenvalue = %d\n",numEigOne);
if numEigOne==nA
   % fprintf("   the closed-loop networked system is stable\n\n");
elseif numEigOne>nA
   fprintf("      the closed-loop networked system is UNstable\n\n");
   stableFlag=false;
else
    fprintf("      Unusual Laplacian matrix. Please check the network\n\n");
    stableFlag=false;
end
indOut1=abs(eigAnetc)>(1+1e-7);
numEigUnstable=sum(indOut1);
if numEigUnstable > 0
    % the close form of the NCS has unstable eigenvalue
    fprintf("     There are %d eigenvalues that are outside unit circle \n", numEigUnstable);
    warning("     Unstable eigenvalue of the networked closed loop system \n");
    stableFlag=false;
end

if stableFlag==true
    disp("     stable system (all eigenvalue <=1)");
end

if DEBUGPLOT==true
    figure('name', 'eigen values of closed-loop networked system'); 
    if stableFlag==true
        strtitle='Localton of eigenvalues (all are stable)'; 
    else
        strtitle='Localton of eigenvalues (with unstable poles)';
    end
    theta = linspace(0,2*pi,300); 
    plot(cos(theta),sin(theta),'g','LineWidth',2); hold on; grid on;
    plot(eigAnetc,'*b'); title(strtitle);
    plot(eigAnetc(indOne),'or');
    plot(eigAnetc(indOut1),'^r','MarkerSize','5');
    axis equal;
end

end

