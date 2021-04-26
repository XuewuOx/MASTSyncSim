clear all
A=[1 1;0 1];
B=[1 1;0 1];
load Ltree20

Ycoptmin=1;
for kk=1:1000
    % K0=rand(2,2);
    p=(rand(2,1)-0.5);
    K0 = place(A,B,p)
    
    [stableFlag, eigAnetc,numEigOne] = chkEigAc(A,B,K0,L);
    fprintf(" ** kk=%d \n",kk);
    if stableFlag==false
        warning("Non-stable control gain matrix");
        continue;
    end
    Linearsearch20shu;
    
    [stableFlag, eigAnetc,numEigOne] = chkEigAc(A,B,K,L);
    if stableFlag==false
        warning("Non-stable control gain matrix");
         continue;
    end
    
    if (Ycopt(end)<1e-9)
      % nID=input('Do want to save the K0 and the optimmized K?? \n','s');
      if Ycopt(end)<Ycoptmin
        ss=num2str(Ycopt(end));
        filename=['goodK_' ss(1) ss(3:5) 'e_' ss(end-1:end)];
    % if (nID=='y' ||  nID=='Y')
        clear sz; % to avoid overwrite sz in main function 
        save(filename);
        Ycoptmin=Ycopt(end);
    end
end
end