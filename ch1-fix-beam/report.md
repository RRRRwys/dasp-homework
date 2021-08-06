



《Fundamentals of Signal Enhancement and Array Signal Processing》7.1-7.5.1学习笔记

# 固定波束

一个固定波束系统就是一个空间滤波器，它可以目标信号的方向形成一个主波束。并且，一可以将一些没有有效信息的干扰信号方向置0。它不会使用从阵列采集的数据或者目标/噪声信号的统计信息。其结果就是滤波器的参数是固定的，阵列的表现不依赖于外部环境的变换。然而，固定波束会通过导向矢量的形式，使用传感器的位置信息，目标方向和干扰方向的信息。因此，需要知道阵列的几何形式。本章主要研究均匀线阵(ULAs)。


## 信号模型和问题形式

我们考虑一个远场平面波，它距离阵列足够远，在消声环境中传播，声速为c(c = 340 m/s)，在一个由M个全指向麦克风构成的均匀线阵上实验。两个相邻的麦克风距离都是 $\delta$ 目标方向是与水平夹角 $\theta$​​。给定导向矢量：
![image-20210803110543609](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110543609.png)


![image-20210803110610005](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110610005.png)

其中j是虚数单位，f是采样频率，$\tau$ 是两个相邻麦克风在方向角为0的时候的延迟。我们定义$\omega = 2 \pi f$是角频率。因为$\cos \theta$是偶函数，所以$d(f,\cos\theta)$也是偶函数。因此，这就把角度限定在了$\theta \in [0,\pi]$。

假定目标信号来自方向角$\theta_d$，我们可以得到它的信号矢量，观测到的信号矢量（长度为$M$）：


![image-20210803110627192](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110627192.png)

其中$Y_m(f)$是第m个麦克风接收的信号，$x(f) = d(f,\cos\theta)X(f)$，$X(f)$是目标信号，$d(f,\cos\theta_d)$是指向为$\theta_d$的导向矢量 （指向目标源），$v(f)$是加性噪声矢量。这样来定义$y(f)$。接下来$y(f)$的自相关矩阵：


![image-20210803110638543](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110638543.png)

$\phi_X(f)$是关于$X(f)$的变量，$\phi_v(f)$是$v(f)$的自相关矩阵。

我们这一章的目标就是设计一个不依赖于统计信号的波束，可以在目标方向角 $\theta_d$ 形成一个主波束，不失真的提取目标信号，并且抑制其他方向的信号。

## 线阵模型

通常，阵列处理或波束是通过在每个麦克风上加一个延迟滤波器，并将滤波后的信号相加。在频域，这等价于给每个麦克风的输出，加一个复值权重，最后求和：
![image-20210803110652997](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110652997.png)
其中$Z(f)$是波束的输出信号，
 ![image-20210803110710198](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110710198.png)
是波束的权重向量，即在频率f上的空域滤波器，
![image-20210803110721626](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110721626.png)
是滤波器的目标信号
![image-20210803110734013](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110734013.png)
是残余噪声。

由于（7.4）的右边两项是不相关的，因此有：
![image-20210803110746138](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110746138.png)
其中


![image-20210803110800344](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110800344.png)
在本文的固定波束中，我们希望获得一个不失真的约束：
![image-20210803110811892](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110811892.png)
这意味着，任何来自目标方向的信号是不失真的。因此，所有的波束都会将（7.11)作为约束。

## 表现评估

在固定波束中，我们通常只关注窄带信号的评估表现，我们将第一个麦克风的信号作为参考信号。

每个波束都有一个对于方向的灵敏度模式：它可以表示对来自不同方向的声音的灵敏度。波束模式或者方向模式描述了波束对于从角度$\theta$入射的平面波，作用在阵列上时的灵敏度。形式化的，它定义为：
![image-20210803110827546](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110827546.png)
通常，$|\mathcal{B}[h(f),\cos\theta]|^2$是功率模式，使用极坐标系绘图。

窄带输入SNR是：
![image-20210803110838833](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110838833.png)
其中$\phi_{V_1}(f) = E[|V_1(f)|^2]$是关于$V_1(f)$的变量，他是$v(f)$的第一个元素。窄带输出SNR定义：
![image-20210803110849124](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110849124.png)
其中

![image-20210803110905986](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110905986.png)

是$v(f)$的互功率谱矩阵。根据之前SNR的定义，我们可以化简得到阵列的增益：
![image-20210803110920780](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110920780.png)

窄带白噪声增益（WNG），是一种最方便的评估阵列缺陷的方案，根据定义（7.16）中的$\Gamma_v(f) = I_M$，其中$I_M$是单位阵：
![image-20210803110931819](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110931819.png)
使用Cauchy-Schwartz不等式：
![image-20210803110940871](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110940871.png)
容易化简得：
![image-20210803110953174](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803110953174.png)
所以，最大WNG是：
![image-20210803111004537](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111004537.png)
他是与频率无关的。令
![image-20210803111015267](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111015267.png)
是2个向量$d(f,\cos\theta)$和$h(f)$之间夹角的余弦值，其中$\Vert·\Vert^2$是$h(f)$​范数。假设不失真约束，我们可以改写WNG为：
![image-20210803111025295](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111025295.png)

另一个重要的指标是麦克风阵列在混响环境下的表现，即窄带指向因子DF，定义为：
![image-20210803111041093](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111041093.png)
其中
![image-20210803111057754](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111057754.png)

我们可以证明，矩阵$\Gamma_{0,\pi}(f)$是：
![image-20210803111106847](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111106847.png)
其中$|\Gamma_{0,\pi}(f)|_{mm} = 1, m = 1,2, \dots, M$。通过Cauchy-Schwartz不等式
![image-20210803111118509](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111118509.png)
可以根据（7.23）发现：
![image-20210803111129252](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111129252.png)
最终，最大指向银子DF就是：
![image-20210803111138597](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111138597.png)
他是关于频率和目标信号夹角的。令
![image-20210803111147828](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111147828.png)
是两个向量夹角的余弦值，假设满足不失真约束，我们可以将DF写作：
![image-20210803111159399](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111159399.png)

## 空间混叠

这里讨论阵列信号处理中遇到的空间混叠问题。它类似于以小于2倍最高频率的采样频率，对时域信号进行采样时，发生的时域混叠现象。
令$\theta_1$​和$\theta_2$​是2分不同的夹角，$\theta_1 \neq \theta_2$​。空间混叠出现在，当$d(f,\cos \theta_1) = d(f, \cos \theta_2)$​​时，即表示声源的方位不明确。令：
![image-20210803111240707](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111240707.png)
或者，等价于
![image-20210803111247677](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111247677.png)
因为
![image-20210803111256867](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111256867.png)

![image-20210803111409405](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111409405.png)

意味着发生空间混叠。

因为$|\cos\theta| \le 1$，所以
![image-20210803111423996](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111423996.png)
等价于
![image-20210803111435851](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111435851.png)
根据（7.32），我们要保证
![image-20210803111445695](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111445695.png)
这就是经典的窄带混叠条件。

## 延迟相加

最著名的固定波束就是延迟相加（delay-and-sum, DS），它可以最大化 WNG：
![image-20210803111455037](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111455037.png)
我们容易得到最优波束滤波器：
![image-20210803111503083](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111503083.png)
使用这个波束时，WNG 和 DF 分别为：
![image-20210803111513399](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111513399.png)
和
![image-20210803111522964](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111522964.png)
因此有
![image-20210803111539078](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111539078.png)
于是有$DF \ge 1$​。当DS波束最大化WNG时，它不会增大扩散场噪声。

我们有波束模式：
![image-20210803111554409](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111554409.png)
DS的波束模式与频率强相关。

另一个表示（7.41）的方式是
![image-20210803111603432](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111603432.png)
其中
![image-20210803111613400](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111613400.png)
它是向量$\Gamma_{0,\pi}^{1/2}(f)d(f,\cos\theta)$和$\Gamma_{0,\pi}^{-1/2}(f)d(f,\cos\theta)$，令$\sigma_1(f)$和$\sigma_M(f)$是$\Gamma_{0,\pi}(f)$​最大和最小的特征值。使用 Kantorovich 不等式：
![image-20210803111628754](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111628754.png)
我们化简得
![image-20210803111640209](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111640209.png)

例子：考虑一个均匀线阵，M个麦克风，如图7.1。假定目标信号入射方向为$\theta$。图7.2 使用了不同的麦克风数量M，绘制了WNG关于频率的变化图。图7.3 展示了不同数量的麦克风，和几种指向角度以及间距。随着麦克风数量增加，WNG和DF都在增加，图7.4 - 7.6展示了M=8时，几种情况的波束模式。主波束指向目标方向。随着频率的增加，主波束的宽度在减小。当$\delta / \lambda$增加，我们可能观察到空间混叠的现象，如图 7.6d

![image-20210803111715165](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111715165.png)

![image-20210803111845197](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111845197.png)

![image-20210803111922314](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111922314.png)

![image-20210803111946087](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803111946087.png)

![image-20210803112007291](C:\Users\YanSongWang\AppData\Roaming\Typora\typora-user-images\image-20210803112007291.png)
