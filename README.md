# Matlab-java
Filtfilt and Spline in java

Filtfilt和Spline函数与Matlab校正后信息一样。
Pyulear函数与Matlab结果不一样，但是函数形态相似。

PyulearBurgC文件的说明-用Burg算法计算信号的功率谱。
1.应用：
PyulearBurgC pyulear = new PyulearBurgC(data, len);
double[] redata = PyulearBurgC.freqMe(data, len, true);
参数：
data：计算功率谱的信号
len:输出数据长度

2.函数是把信号的功率谱投射在0-0.5的区间上。
因此实际输出信号对应的频率f = n/len*fs/2;
fs为信号采样率。n为从1：len的遍历。
