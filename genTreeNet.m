function [NetTree]= genNet(L)

%     clear all;
%     close all;
%     clc;
%     load('OriginalDirectedGraphwith50Nodes&60Edges.mat','L');

%%
% node 23 -> node 1 
L(23,23) = 1;
for i = 24:50
    L(23,i) = 0;
end
sum(L(23,:));

% node 45 -> node 23 
L(45,45) = 1;
L(45,7) = 0;
sum(L(45,:));

%%
% node 35 -> node 23 
L(35,35) = 1;
L(35,3) = 0;
L(35,12) = 0;
L(35,36) = 0;
sum(L(35,:));

% node 12 -> node 35 
L(12,12) = 1;
L(12,10) = 0;
sum(L(12,:));

% node 10 -> node 12 
L(10,10) = 1;
L(10,13) = 0;
L(10,19) = 0;
sum(L(10,:));

% node 13 -> node 10 
L(13,13) = 1;
L(13,21) = 0;
sum(L(10,:));

% node 21 -> node 13 
L(21,21) = 1;
L(21,4) = 0;
sum(L(21,:));

% node 4 -> node 21 
L(4,4) = 1;
L(4,18) = 0;
L(4,32) = 0;
sum(L(4,:));

% node 32 -> node 4 
L(32,32) = 1;
L(32,17) = 0;
sum(L(32,:));

% node 17 -> node 32 
L(17,17) = 1;
L(17,27) = 0;
sum(L(17,:));

% node 18 -> node 4 
L(18,18) = 1;
L(18,31) = 0;
sum(L(18,:));

%%
% node 3 -> node 35 
L(3,3) = 1;
L(3,5) = 0;
L(3,50) = 0;
sum(L(3,:));

% node 50 -> node 3 
L(50,50) = 1;
L(50,38) = 0;
L(50,34) = 0;
sum(L(50,:));

% node 34 -> node 50 
L(34,34) = 1;
L(34,37) = 0;
sum(L(34,:));

% node 37 -> node 34 
L(37,37) = 1;
L(37,28) = 0;
L(37,44) = 0;
L(37,46) = 0;
sum(L(37,:));

% node 28 -> node 37 
L(28,28) = 1;
L(28,20) = 0;
sum(L(28,:));

% node 20 -> node 28 
L(20,20) = 1;
L(20,2) = 0;
L(20,42) = 0;
sum(L(20,:));

% node 42 -> node 20
L(42,42) = 1;
L(42,6) = 0;
sum(L(42,:));

%%
% node 5 -> node 3
L(5,5) = 1;
L(5,11) = 0;
sum(L(5,:));

% node 11 -> node 5
L(11,11) = 1;
L(11,8) = 0;
L(11,14) = 0;
L(11,33) = 0;
sum(L(11,:));

% node 14 -> node 11
L(14,14) = 1;
L(14,43) = 0;
sum(L(14,:));

% node 8 -> node 11
L(8,8) = 1;
L(8,29) = 0;
L(8,49) = 0;
sum(L(8,:));

% node 29 -> node 16
L(29,29) = 1;
L(29,16) = 0;
sum(L(29,:));

% node 33 -> node 11
L(33,33) = 1;
L(33,9) = 0;
sum(L(33,:));

% node 9 -> node 33
L(9,9) = 1;
L(9,15) = 0;
L(9,24) = 0;
L(9,30) = 0;
sum(L(9,:));

% node 24 -> node 9
L(24,24) = 1;
L(24,25) = 0;
sum(L(24,:));

% node 25 -> node 24
L(25,25) = 1;
L(25,22) = 0;
sum(L(25,:));

% node 22 -> node 25
L(22,22) = 1;
L(22,41) = 0;
L(22,48) = 0;
sum(L(22,:));

% node 41 -> node 22
L(41,41) = 1;
L(41,40) = 0;
sum(L(41,:));

% node 40 -> node 41
L(40,40) = 1;
L(40,39) = 0;
sum(L(40,:));

NetTree = L;

%     % directed graph
%     d=diag(NetTree);  % d is the dialog vector of matrix L 
%     A=diag(d)-L;
%     netG=digraph(A,'omitselfloops');
%     figure('name', 'Network Topology'); plot(netG,'LineWidth',2);
%     title("a directed  graph ");

end