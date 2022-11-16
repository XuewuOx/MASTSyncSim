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

2021.04.22
完成基于LMI 线性搜索的增益矩阵K的设计，采用了第一种模式， designK() 由脚本 Linearsearch20shu.m 
及三个setLMIProb_****.m 来完成，返回优化后的增益矩阵K. 目前仅针对特定的20节点的树形拓扑设计

## Multi-agent Systems for Time Synchronisation in Wireless Networks

- Whenever using any part of the code in the branch [Zong2021c_TII](https://github.com/zongyan/MASTSyncSim/tree/Zong2021c_TII), please cite the following [paper](http://yzong.com/pub/Zong2022b.pdf)

  - Y. Zong, S. Liu, X. Liu, S. Gao, X. Dai, and Z. Gao., “Robust Synchronised Data Acquisition for Biometric Authentication,” IEEE Trans. Ind. Informat., DOI: 10.1109/TII.2022.3182326

- Utilising any part of the code in the branch [Zong2021b_IoTJ](https://github.com/zongyan/MASTSyncSim/tree/Zong2021b_IoTJ), please cite the following [paper](http://yzong.com/pub/Zong2022a.pdf)

  - Y. Zong, X. Dai, Z. Wei, M. Zou, W. Guo, and Z. Gao., “Robust Time Synchronisation for Industrial Internet of Things by $H_\infty$ Output Feedback Control,” IEEE Internet of Things J., DOI: 10.1109/JIOT.2022.3144199

- If Using any part of the code in the branch [Zong2021a_IoTJ](https://github.com/zongyan/MASTSyncSim/tree/Zong2021a_IoTJ), please cite the following [paper](http://yzong.com/pub/Zong2021.pdf)

  - Y. Zong, X. Dai, S. Gao, P. Canyelles-Pericas, and S. Liu., “PkCOs: Synchronisation of Packet-Coupled Oscillators in Blast Wave Monitoring Networks,” IEEE Internet of Things J., vol. 9, no. 13, pp. 10862-10871, Jul. 2022 DOI: 10.1109/JIOT.2021.3126059
