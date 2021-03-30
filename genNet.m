function [netG netTree]= genNet(nNode,nEdge,isCmdPromtEnabled)
% genNet 
% to generate a random network 
% Input arguments: nNode   expected number of nodes, default 5 nodes
%                  nEdge   expected number of edges, default 15 edges
%                  isCmdPromtEnabled   1 allow user select the suitable network
%                  graph
% Output arguments: netG   the graph object of the network generated
            
% if command prompt is enabled 
 if nargin==3
     cmdp= isCmdPromtEnabled==1;
 else
     cmdp=false
 end
 if nargin==1
        if nNode>4
            nEdge=2*nNode;
        else
            nEdge=6;
        end
 elseif nargin==0
     nNode=5; nEdge=8;
 end
 
fh=figure('Name', 'Network Graph');
while true  % for user decide if a graph is suitable 
    for i=50:-1:1   % search for a connected graph
        % generate nEdge connections randommanly
        st = randi(nNode, nEdge, 2);
        
        % check and remove selfloop and duplicated edges
        nSelfloop=1;
        nDupEdge=1;
        while (nSelfloop ~= 0 || nDupEdge~=0 )
            % check self loop
            selfLp=diff(st,[],2)==0;
            nSelfloop=sum(selfLp);
            if nSelfloop ~=0
                st(selfLp, :)=[]; % remove self loop
            end
            % to check duplicated edge
            % sort the edges according to the sum and different of the edge's node IDs
            st2=[sum(st,2), abs(diff(st,[],2)), st];
            st2=sortrows(st2);
            st3=diff(st2(:,1:2),[],1);
            % check if the edge is the same as the previous one
            % if the sum and abs(diff) of the edge's two node IDs are the same
            %  then this is a duplicated edge. remove
            dupEdge=and(st3(:,1)==0, st3(:,2)==0);
            dupEdge=[false;dupEdge];
            nDupEdge=sum(dupEdge);
            st2(dupEdge,:)=[]; % remove repeated edge
            % restore st
            st=st2(:,3:4);
            
            %replace these removed edges
            st = [st;randi(nNode,nSelfloop+nDupEdge,2)];
            fprintf("found and replaced %d selfloops, %d duplicated Edges\n",nSelfloop,nDupEdge);
        end
        
       
        
        fprintf("No selfloops, neither duplicated edges\n");
        s = st(:,1); t = st(:,2);
        netG=graph(s,t);
        if ismultigraph(netG)  % double check duplicated edges
            warning("the created graph has multiple duplicated edges between any two nodes. Please check your network graph generation codes");
        end
        L=laplacian(netG);
        
        % check if is L is connected
        % Get eigenvalues of laplacian
        lambda = eig(L);
        % Get number of zero eigenvalues
        num=sum(lambda<=0.00001);
        discon=num>1; % more than one 0 eigenvalue, disconnected
        if discon
            disp (['Graph is not connected. ' num2str(num) ' connected subgraphs.']);
            fprintf('+++++++   Try again (i=%d) +++++++\n',i);
        else
            disp('Graph is connected!');
            break;
        end
    end % end of for loop, search connected graph

    figure(fh); hold off;
    hfig=plot(netG);title(sprintf("network of %d nodes, %d edges",nNode, nEdge));
    
    [netTree,pred] = minspantree(netG);
    highlight(hfig,netTree)

    fprintf("the Laplacian Matrix of the generated network is \n")
    L=laplacian(netG); % a sparse matrix of L
    % disp(L); % diaplay spares matrix L
    disp(full(L)); % convert sparse matrix to full storage

    if i<=1
        warning("Failed to generate a connected network with %d nodes and %d edges",...
            nNode, nEdge);
        disp("Press any key to continue ...");
        pause
    else
        fprintf("A connected network of %d node, %d edges are created successfully\n",nNode, nEdge);
    end
    if ~cmdp
        break;
    end
    
    disp('Is the auto generated network ok (Y/R/Q)? \n');
    disp('    yY  ---   grah is good. Process to next step');
    disp('    rR  ---   retry to generate another graph');
    disp('    qQ  ---   quit without network graph. return with []');
    x = input( 'Type your choice:', 's');
    if x == 'y' || x=='Y'
        fprintf("Both graph (%d nodes, %d edges) and its minimal spanning tree (%d edges) are created\n", ...
            nNode,nEdge,nEdge);
        break;
    elseif x=='r' || x=='R'
        disp('Try again to generate another network...\n');
        % close(fh);fh=[];
        continue;
    else
        netG=[];
        disp('NO graph is created. User select to quit.\n');
        return;
    end
end % end while
    

end

