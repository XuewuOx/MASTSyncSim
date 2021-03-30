# MASTSyncSim
Time Synchronization by Multi-agent System Simulation

The initial version has the following features:
(1) Create an arbitrary graph or one from a predefined Laplacian matrix. 
(2) Allow the user to decide using a graph or tree for Time Sync
(3) Allow the user to set reference clock(s) or non-servo clocks
(4) Simulaiton animation and results visualization.  

控制器的m函数文件有两种调用方案：
（1） 比例控制器增益矩阵设计函数 designK（）： 用于确定控制器增益矩阵， 调用designK（）时，返回控制增益矩阵K
(2) 控制器函数 getCtrlISig(x): 输入参数为网络当前状态向量，输出参数为 控制输入信号的值u
