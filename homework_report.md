# 均衡器

# webrtc aec

核心代码：https://android.googlesource.com/platform/external/webrtc/+/e48d5845c8b35de2ab73ea055c18a61fa3a9f0be/src/modules/audio_processing/aec/

新版webrtc中已经不再使用，目前是 AEC3

主要包括延迟估计、自适应滤波、非线性处理三个部分。

提供了matlab代码 fullaec.m，其中不含延迟估计部分

## 算法介绍

最速梯度下降
h(n+1) = h(n) - g(n) miu(n)
g(n) 是均方误差对于h(n)的导数，miu(n)是步长
LMS
使用瞬时均方误差，代替均方误差，可得：
(d - x * h)^2 = d^2 - 2 * d * x * h + x * h * x * h
g1(n) = - 2 * d * x + 2 * x^2 * h = - 2 (d - x * h) * x = - 2 * e * x
可得
h(n+1) = h(n) + \miu e(n) x(n)

NLMS
x(n)较大时，会引起梯度放大，较小时收敛较慢，因此 \miu 对瞬时能量归一化。
h(n+1) = h(n) + miu(n) e(n) x(n)
miu(n) = \mu / |x(n)|^2

NLMS 对于较大的延时，需要抽头多，延时50ms以上时，按照16kHz语音采样率来看，抽头 > 16000 * 0.05(800) 是最低的要求，如果希望覆盖反射回声，那滤波器的抽头就得数以千计了。

PBFDAF
...

思想上与NLMS基本一致，不同点在于用之前的多个块，一起对当前时刻进行估计

## 整体框架介绍

## 代码分析

~~通过阅读c代码了解延迟估计部分~~

通过分析fullaec.m了解自适应滤波和非线性处理部分

PBFDAF 公式

## 实验结果


