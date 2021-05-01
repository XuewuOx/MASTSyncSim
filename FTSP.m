% % ***********************************************************************
% * Descriptoin: this file is used to simulate the performance of FTSP in 
%               a network, and the least squares regression is used in this
%               work
% % ***********************************************************************
%
% comments: 1. the buffer size of linear regression in FTSP is 8.
%           2. only synchronsing with the root node or synchronised node
%           3. the buffer of linear regression is only cleared when K number 
%              of synchronisation periods, where K is pre-definied. This 
%              clearing buffer feature can be ignored now.

% the pseudo-code of FTSP is as follows
%
% node i receives wireless packets from node j
% if node j is the root node or synchronised node
%   node i records the clock offset
% else 
%   node i diacards the received clock offset
% end if
% 
% if consistent eight offsets are obtained 
%   estimating the i-th skew via linear regression, and also i-th offset 
%   corrected the skew and offset of node i
% end if



