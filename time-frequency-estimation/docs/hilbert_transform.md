# Hilbert 变换

## 1. 方法目标

Hilbert 变换常用于构造解析信号。对于单分量窄带实信号，解析信号可以分离出包络和相位，从而进一步估计瞬时频率。

普通实信号 $x(t)$ 只有一个实部。通过 Hilbert 变换可以构造虚部：

$$
z(t)=x(t)+j\mathcal{H}\{x(t)\}
$$

其中 $z(t)$ 称为解析信号。

## 2. 数学原理

Hilbert 变换定义为：

$$
\mathcal{H}\{x(t)\}=\frac{1}{\pi}\,\mathrm{p.v.}\int_{-\infty}^{+\infty}\frac{x(\tau)}{t-\tau}\,d\tau
$$

其中 $\mathrm{p.v.}$ 表示柯西主值积分。

在频域中，Hilbert 变换等价于对正负频率施加不同相位偏移：

$$
\mathcal{F}\{\mathcal{H}[x(t)]\}=-j\,\mathrm{sgn}(f)X(f)
$$

解析信号为：

$$
z(t)=x(t)+j\mathcal{H}\{x(t)\}
$$

解析信号可以写成极坐标形式：

$$
z(t)=A(t)e^{j\phi(t)}
$$

其中：

$$
A(t)=|z(t)|
$$

是瞬时包络，

$$
\phi(t)=\arg z(t)
$$

是瞬时相位。

对于幅度调制窄带信号：

$$
x(t)=A(t)\cos(2\pi f_ct)
$$

Hilbert 解析信号的模长可以较好恢复包络 $A(t)$。

## 3. 与代码的对应关系

示例脚本：`examples/hilbert_transform_example.py`

`make_am_signal(fs, duration)` 构造幅度调制信号：

$$
x(t)=A(t)\cos(2\pi f_ct)
$$

其中：

$$
A(t)=1+0.45\sin(2\pi\cdot 3t)
$$

载波频率为：

$$
f_c=80\ \mathrm{Hz}
$$

`analytic_signal(x)` 调用：

```python
z = signal.hilbert(x)
```

对应：

$$
z(t)=x(t)+j\mathcal{H}\{x(t)\}
$$

代码中的：

```python
envelope = np.abs(z)
```

对应：

$$
A(t)=|z(t)|
$$

代码中的：

```python
phase = np.unwrap(np.angle(z))
```

对应先计算：

$$
\phi(t)=\arg z(t)
$$

再进行相位展开。相位展开用于去除 $\arg$ 函数导致的 $2\pi$ 跳变，使相位随时间连续增长。

第三幅图绘制：

```python
np.real(z), np.imag(z)
```

它展示解析信号在复平面上的轨迹，实部是原始信号，虚部是 Hilbert 变换结果。

## 4. 运行方式

```bash
python3 examples/hilbert_transform_example.py
```

输出图片：

```text
outputs/hilbert_transform_example.png
```
