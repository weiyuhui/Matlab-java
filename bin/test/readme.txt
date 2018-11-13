Spline_Mat.java：实现了 Matlab的spline函数，经过比较结果完全一致；
testJama.java：是测试用例。

说明：
由于Matlab的spline函数用到了稀疏矩阵，因此引用了ujmp-complete-0.3.0.jar这个库。
说明文档及下载详情见：https://ujmp.org/

缺点：
由于ujmp-complete-0.3.0.jar库不支持维度不同的矩阵除法，而spline函数的实现依赖矩阵的除法，因而本文采取的方法是先求逆矩阵，然后逐行求乘法再相加，有可能遇到超大矩阵会内存溢出。
