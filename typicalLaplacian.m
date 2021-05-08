% typical laplacian matrices 
   L3=[1 -1 0;
        -1 2 -1;
        0 -1 1];
    L4=[1 -1 0 0;
        -1 2 -1 0;
        0 -1 2 -1;
        0 0  -1 1];
    
    L4r=[1 -1 0 0;
        -1 3 -1 -1;
        0 -1 1 0;
        0 -1 0 1];
    
L1way= ...
  [  1	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	0	0	2	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	0	0	0	2	0	0	0	0	0	-1	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	1	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	1	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	1	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	1	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	1	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	1	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	1	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	1];


L2=[  6	-1	-1	-1	-1	-1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	0	0	7	-1	-1	-1	-1	-1	-1	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	0	0	0	0	0	7	-1	-1	-1	-1	-1	-1	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	1	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	1	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	1	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	6	-1	-1	-1	-1	-1;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	1	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	1	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	1	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	1];


 Lshu20=[2 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
      -1 3 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      -1 0 3 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1 0 3 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0;
      0 -1 0 0 3 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0;
      0 0 -1 0 0 3 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0;%row 6?
      0 0 -1 0 0 0 3 0 0 0 0 0 0 -1 -1 0 0 0 0 0;
      0 0 0 -1 0 0 0 3 0 0 0 0 0 0 0 -1 -1 0 0 0;
      0 0 0 -1 0 0 0 0 3 0 0 0 0 0 0 0 0 -1 -1 0;
      0 0 0 0 -1 0 0 0 0 2 0 0 0 0 0 0 0 0 0 -1;% row?10
      0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 -1 -0 0 0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0 0 0 0 0;%row 15
      0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0 0;
      0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0;
      0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1;
    ];% 树形20个节点双向结构?

Lshu21=[
    2  -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
   -1   2  0  -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
   -1   0  3  0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0   -1 0 3 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0  0 -1 0 3 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0;
    0  0 -1 0 0 3 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0;
    0  0 0 -1 0 0 3 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0;%row 6?
    0  0 0 -1 0 0 0 3 0 0 0 0 0 0 -1 -1 0 0 0 0 0;
    0  0 0 0 -1 0 0 0 3 0 0 0 0 0 0 0 -1 -1 0 0 0;
    0  0 0 0 -1 0 0 0 0 3 0 0 0 0 0 0 0 0 -1 -1 0;
    0  0 0 0 0 -1 0 0 0 0 2 0 0 0 0 0 0 0 0 0 -1;% row?10
    0  0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
    0  0 0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
    0  0 0 0 0 0 -1 -0 0 0 0 0 0 1 0 0 0 0 0 0 0;
    0  0 0 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0;
    0  0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0 0 0 0 0;%row 15
    0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0 0 0 0;
    0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0 0;
    0  0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0;
    0  0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1 0;
    0  0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1;
    ];% 树形21个节点双向结构?
 
% 21-node tree, same as Lshu21, but Node 1 (root) fails
Lshu21_n1F=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,2,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,-1,3,0,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
   0   -1 0 3 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0  0 -1 0 3 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0;
    0  0 -1 0 0 3 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0;
    0  0 0 -1 0 0 3 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0;%row 6?
    0  0 0 -1 0 0 0 3 0 0 0 0 0 0 -1 -1 0 0 0 0 0;
    0  0 0 0 -1 0 0 0 3 0 0 0 0 0 0 0 -1 -1 0 0 0;
    0  0 0 0 -1 0 0 0 0 3 0 0 0 0 0 0 0 0 -1 -1 0;
    0  0 0 0 0 -1 0 0 0 0 2 0 0 0 0 0 0 0 0 0 -1;% row?10
    0  0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
    0  0 0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
    0  0 0 0 0 0 -1 -0 0 0 0 0 0 1 0 0 0 0 0 0 0;
    0  0 0 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 0;
    0  0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0 0 0 0 0;%row 15
    0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0 0 0 0;
    0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0 0;
    0  0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0;
    0  0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1 0;
    0  0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 1;
    ];% 21-node tree, but root node (node 1) failure

   Lshu12=[2 -1 -1 0 0 0 0 0 0 0 0 0;
      -1 3 0 -1 -1 0 0 0 0 0 0 0;
      -1 0 3 0 0 -1 -1 0 0 0 0 0;
      0 -1 0 3 0 0 0 -1 -1 0 0 0;
      0 -1 0 0 3 0 0 0 0 -1 -1 0;
      0 0 -1 0 0 1 0 0 0 0 0 0;
      0 0 -1 0 0 0 2 0 0 0 0 -1;
      0 0 0 -1 0 0 0 1 0 0 0 0;
      0 0 0 -1 0 0 0 0 1 0 0 0;
      0 0 0 0 -1 0 0 0 0 1 0 0;
      0 0 0 0 -1 0 0 0 0 0 1 0;
      0 0 0 0 0 0 -1 0 0 0 0 1
    ];%树形12个节点双向结构?
Lcu20= ...
  [  1	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    -1	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;%绗簲琛?
    -1	0	0	0	0	2	0	0	0	0	-1	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	0	1	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	0	0	0	1	0	0	0	0	0	0	0	0	0	0;%绗崄琛?
    0	0	0	0	0	-1	0	0	0	0	2	0	0	0	0	-1	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	1	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	0	1	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	0	0	0	1	0	0	0	0	0;%绗崄浜旇
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	-1	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	1	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	1	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	1];%绗簩鍗佽

Lhuan20=...
  [  2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1;
    -1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1	0;
    0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2	-1;
    -1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	2;
   ];


    % L12: 12 nodes linear network
    L12=[1 -1 0 0 0 0 0 0 0 0 0 0;
        -1 2 -1 0 0 0 0 0 0 0 0 0;
        0 -1 2 -1 0 0 0 0 0 0 0 0;
        0 0 -1 2 -1 0 0 0 0 0 0 0;
        0 0 0 -1 2 -1 0 0 0 0 0 0;
        0 0 0 0 -1 2 -1 0 0 0 0 0;
        0 0 0 0 0 -1 2 -1 0 0 0 0;
        0 0 0 0 0 0 -1 2 -1 0 0 0;
        0 0 0 0 0 0 0 -1 2 -1 0 0;
        0 0 0 0 0 0 0 0 -1 2 -1 0;
        0 0 0 0 0 0 0 0 0 -1 2 -1;
        0 0 0 0 0 0 0 0 0 0 -1 1;
        ];%单线
    % L12c: 12 nodes circular network
    L12c=[2 -1 0 0 0 0 0 0 0 0 0 -1;
        -1 2 -1 0 0 0 0 0 0 0 0 0;
        0 -1 2 -1 0 0 0 0 0 0 0 0;
        0 0 -1 2 -1 0 0 0 0 0 0 0;
        0 0 0 -1 2 -1 0 0 0 0 0 0;
        0 0 0 0 -1 2 -1 0 0 0 0 0;
        0 0 0 0 0 -1 2 -1 0 0 0 0;
        0 0 0 0 0 0 -1 2 -1 0 0 0;
        0 0 0 0 0 0 0 -1 2 -1 0 0;
        0 0 0 0 0 0 0 0 -1 2 -1 0;
        0 0 0 0 0 0 0 0 0 -1 2 -1;
        -1 0 0 0 0 0 0 0 0 0 -1 2;
        ];%环形
    % Lsu12=[2 -1 -1 0 0 0 0 0 0 0 0 0;
    %   -1 3 0 -1 -1 0 0 0 0 0 0 0;
    %   -1 0 3 0 0 -1 -1 0 0 0 0 0;
    %   0 -1 0 3 0 0 0 -1 -1 0 0 0;
    %   0 -1 0 0 3 0 0 0 0 -1 -1 0;
    %   0 0 -1 0 0 1 0 0 0 0 0 0;
    %   0 0 -1 0 0 0 2 0 0 0 0 -1;
    %   0 0 0 -1 0 0 0 1 0 0 0 0;
    %   0 0 0 -1 0 0 0 0 1 0 0 0;
    %   0 0 0 0 -1 0 0 0 0 1 0 0;
    %   0 0 0 0 -1 0 0 0 0 0 1 0;
    %   0 0 0 0 0 0 -1 0 0 0 0 1
    % ];%树形
