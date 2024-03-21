# Coupling-Matrix-Extraction-Method-With-High-Order-Transmission-Coefficient-Fitting-Strategy
These MATLAB source codes are related to the paper "Application of Higher-Order S21 Fitting Strategy in Filter Coupling Matrix Extraction".

本文件为论文《高阶S21拟合策略在耦合矩阵提取方法中的运用》的配套源代码以及测试用数据。

============================================================================================

论文地址：https://kns.cnki.net/kcms2/article/abstract?v=smPsKIJgVaD7nNc1B3Fr41t1dEO5s5ImDrpCdTjnJcbH4ggR6ttWUyoUuLDclWWNNhjAQ4fuKwIPRUZBAY0FeN3_EufibWq3Xws91gndAXDHLv0Z-bZCB422-PiiU-ePUvdUlptWprI=&uniplatform=NZKPT&language=CHS
##注意

所有代码采用MATLAB语言编写，使用MATLAB R2018b软件进行试验。

脚本中所用的部分MATLAB子函数来源于Communications Toolbox工具箱，请使用版本高于R2018b的、含有Communications Toolbox工具箱的MATLAB 软件运行主函数脚本文件。

============================================================================================
##脚本文件说明

子函数脚本文件说明：

analyseCM.m：分析耦合矩阵S参数

De_Embedding_2p_onlyarg.m：双端口滤波器网络去嵌入

Reads2p.m：读取s2p格式文件

StoY.m：端口阻抗1欧姆下，S参数转换Y参数

to_foldedCM.m：横向耦合矩阵旋转至规范折叠型耦合矩阵

主函数脚本文件说明：

OriginalCauchyMethod.m：传统柯西方法，参考文献：
[7]ZHAO Ping, WU Keli. A New Computer-Aided Tuning Scheme for General Lossy Coupled-Resonator Bandpass Filters Based on the Cauchy Method [J]. HKIE Transactions, 2016, 23(1): 52-62.

CauchyMethod.m：应用高阶S21拟合策略的柯西方法（具体操作步骤可参考《高阶S21拟合策略在耦合矩阵提取方法中的运用》）

MVF.m：基于模型的矢量拟合方法，参考文献：
[11]ZHAO Ping, WU Keli. Model-Based Vector-Fitting Method for Circuit Model Extraction of Coupled-Resonator Diplexers [J]. IEEE Transaction on Microwave Theory and Techniques, 2016, 64(6): 1787-1797.

VFSp.m：通过矢量拟合定位S参数零点的耦合矩阵提取方法（具体操作步骤可参考《高阶S21拟合策略在耦合矩阵提取方法中的运用》）

============================================================================================
##测试数据文件说明

Sx.s2p：通带1920-1980MHz，9阶-3零点同轴滤波器六组状态下的HFSS仿真数据。

滤波器参数：

N = 9; % Order of the filter

Nz = 3; % Number of the transmission zeros

CF = sqrt(1920*1980)*1e6; % Hz center mapping frequency

BW = (1980-1920)*1e6; % Hz mapping bandwidth

1900mea.s2p：通带1920-1980MHz，6阶-2零点悬置带线滤波器测试数据。

滤波器参数：

N = 6; % Order of the filter

Nz = 2; % Number of the transmission zeros 

CF = sqrt(1920*1980)*1e6; % Hz center mapping frequency

BW = (1980-1920)*1e6; % Hz mapping bandwidth

============================================================================================

By YellowBook

Time: 20240107
