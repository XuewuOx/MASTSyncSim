function [netG netTree]= genNet(nNode,nEdge,isCmdPromtEnabled)
% genNet is to generate a random network 
% Input arguments: nNode   expected number of nodes, default 5 nodes
%                  nEdge   expected number of edges, default 15 edges
%                  isCmdPromtEnabled   1 allow user select the suitable network
%                  graph
% Output arguments: netG   the graph object of the network generated
% 
%% 
% Description: according to the value of nargin (which returns number of 
%              function input arguments), determining 'user select the 
%              suitable network', or using the default configurations

% if command prompt is enabled 
 if nargin==3 % the number of function input arguments is three
     cmdp= isCmdPromtEnabled==1;
 else
     cmdp=false
 end
 if nargin==1 % only one input
        if nNode>4
            nEdge=2*nNode;
        else
            nEdge=6;
        end
 elseif nargin==0 % by default, 5 nodes and 8 edges
     nNode=5; nEdge=8;
 end
%% 
% using randi to generate a nEdge-by-2 matrxi; using the differences between 
% columns to determine are there any self-loops; and using the differences 
% between rows, in order to determine are there any self-loops duplicated edges

fh=figure('Name', 'Network Graph');
while true  % for user decide if a graph is suitable 
    for i=200:-1:1   % search for a connected graph
        % generate nEdge connections randommanly
        st = randi(nNode, nEdge, 2); % returns an nEdge-by-2 matrix of 
                                     % pseudorandom integers drawn from the 
                                     % discrete uniform distribution on the 
                                     % interval [1, nNode].
                                     % 这个矩阵，会放到graph，see line 88
        
        % checking and remove self-loop and duplicated edges
        nSelfloop=1;
        nDupEdge=1;
        while (nSelfloop ~= 0 || nDupEdge~=0 )
            % checking self-loops
             
            selfLp=diff(st,[],2)==0;% diff(st,[],2): created a (st.row)-by-2 
                                    % matrix, then compute the first-order 
                                    % ([] represents []-th order) difference 
                                    % between the columns.
                                    
                                    % determining the element is zero not
                                    % not, if yes return logic 'one'

            nSelfloop=sum(selfLp);  % the number of 'one', one represents 
                                    % self-loop, 也代表st中st(i,j)&st(i,j+1)
                                    % 相等的个数
                                    
            if nSelfloop ~=0
                st(selfLp, :)=[]; % remove self-loops, 即删除st中st(i,j)&st(i,j+1)
                                  % 相等的行
            end
            
            % checking duplicated edges
            % sort the edges according to the sum and different of the edge's node IDs
            st2=[sum(st,2), abs(diff(st,[],2)), st];
            st2=sortrows(st2); % sortrows(A) sorts the rows of a matrix in 
                               % ascending order based on the elements in 
                               % the first column. When the first column 
                               % contains repeated elements, sortrows sorts 
                               % according to the values in the next column 
                               % and repeats this behavior for succeeding 
                               % equal values. 整行重新排序
            st3=diff(st2(:,1:2),[],1); % 第一列是“sum和”，第二列是“差的绝对值”
                                       % 第二行减去第一行，依次逐行相减
            % checking if the edge is the same as the previous one. if the 
            % sum and abs(diff) of the edge's two node IDs are the same
            % then this is a duplicated edge. remove
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
        
        fprintf("No self-loops, neither duplicated edges\n");
        s = st(:,1); t = st(:,2);
        netG=graph(s,t); % constructing an undirected graph with edges 
                         % specified by the node pairs (s, t). 
                         
        if ismultigraph(netG)  % determining whether a graph has multiple edges. 
                               % double check duplicated edges
            warning("the created graph has multiple duplicated edges between any two nodes. Please check your network graph generation codes");
        end
        
        %% 
        % according the eigenvalue of a graph's laplacian matrix, to determine
        % are all nodes in a graph connected
        % 
        % ToDo: what is the difference between enginvalues and laplacian in
        % a graph
        
        L=laplacian(netG);        
        
        % checking if is L is connected
        lambda = eig(L); % Get eigenvalues of laplacian
        num=sum(lambda<=0.00001); % Get number of zero eigenvalues
        discon=num>1; % more than one 0 eigenvalue, disconnected
        if discon
            disp (['Graph is not connected. ' num2str(num) ' connected subgraphs.']);
            fprintf('+++++++   Try again (i=%d) +++++++\n',i);
        else
            disp('Graph is connected!');
            break;
        end
    end % end of for loop, search connected graph
%% generating an undirected graph (could be minimum spanning tree, could be graph)
    figure(fh); hold off;
    hfig=plot(netG);title(sprintf("network of %d nodes, %d edges",nNode, nEdge));
    
    [netTree,pred] = minspantree(netG); % minimum spanning tree of undirected graph
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
        fprintf("Both undirected graph (%d nodes, %d edges) and its minimal spanning tree (%d edges) are created\n", ...
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

