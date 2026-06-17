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

这里必须使用柯西主值，是因为当 $\tau=t$ 时，积分核：

$$
\frac{1}{t-\tau}
$$

有奇点。普通反常积分会把奇点左右两侧分开考察，但 Hilbert 变换需要保留这个奇异核的“奇对称抵消”性质。因此更严格地写：

$$
\mathcal{H}\{x(t)\}
=
\frac{1}{\pi}
\lim_{\varepsilon\to 0^+}
\left[
\int_{-\infty}^{t-\varepsilon}
\frac{x(\tau)}{t-\tau}\,d\tau
+
\int_{t+\varepsilon}^{+\infty}
\frac{x(\tau)}{t-\tau}\,d\tau
\right]
$$

也可以令 $u=t-\tau$，写成卷积形式：

$$
\mathcal{H}\{x(t)\}
=
x(t) * \frac{1}{\pi t}
$$

但这里的 $\frac{1}{t}$ 不是普通函数意义下的卷积核，而是按柯西主值理解的广义函数：

$$
\frac{1}{\pi t}
\quad\Longleftrightarrow\quad
\frac{1}{\pi}\,\mathrm{p.v.}\frac{1}{t}
$$

在频域中，Hilbert 变换等价于对正负频率施加不同相位偏移：

$$
\mathcal{F}\{\mathcal{H}[x(t)]\}=-j\,\mathrm{sgn}(f)X(f)
$$

也就是说：

$$
\mathcal{H}\{x(t)\}
\quad\Longleftrightarrow\quad
\begin{cases}
-jX(f), & f>0 \\
0, & f=0 \\
+jX(f), & f<0
\end{cases}
$$

所以 Hilbert 变换不会改变频率幅度，只改变相位：

- 正频率分量相位移动 $-\frac{\pi}{2}$。
- 负频率分量相位移动 $+\frac{\pi}{2}$。
- 直流分量被置为 $0$。

解析信号为：

$$
z(t)=x(t)+j\mathcal{H}\{x(t)\}
$$

令 $Z(f)=\mathcal{F}\{z(t)\}$，则：

$$
\begin{aligned}
Z(f)
&=
X(f)+j\mathcal{F}\{\mathcal{H}[x(t)]\} \\
&=
X(f)+j[-j\,\mathrm{sgn}(f)X(f)] \\
&=
[1+\mathrm{sgn}(f)]X(f)
\end{aligned}
$$

因此：

$$
Z(f)=
\begin{cases}
2X(f), & f>0 \\
X(0), & f=0 \\
0, & f<0
\end{cases}
$$

这就是“解析信号”的关键：它把实信号中的负频率分量去掉，并把正频率分量加倍。实信号本来具有共轭对称频谱，正负频率携带的是同一份信息；解析信号把这份信息集中到正频率一侧，因此可以自然地写成复指数形式。

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

Hilbert 解析信号的模长可以较好恢复包络 $A(t)$。其原因可以从窄带调幅信号推导出来。

先把余弦写成复指数：

$$
x(t)
=
A(t)\cos(2\pi f_ct)
=
\frac{1}{2}A(t)e^{j2\pi f_ct}
+
\frac{1}{2}A(t)e^{-j2\pi f_ct}
$$

如果 $A(t)$ 是慢变化包络，它的频谱主要集中在低频。设 $A(t)$ 的最高主要频率远小于载波频率 $f_c$，则 $A(t)e^{j2\pi f_ct}$ 的频谱主要位于正频率 $+f_c$ 附近，$A(t)e^{-j2\pi f_ct}$ 的频谱主要位于负频率 $-f_c$ 附近，两者基本不重叠。

由上面的解析信号频域性质可知，构造 $z(t)$ 时会去掉负频率、加倍正频率：

$$
z(t)
=
2\cdot \frac{1}{2}A(t)e^{j2\pi f_ct}
=
A(t)e^{j2\pi f_ct}
$$

于是：

$$
|z(t)|
=
|A(t)e^{j2\pi f_ct}|
=
|A(t)|\cdot |e^{j2\pi f_ct}|
=
|A(t)|
$$

通常包络 $A(t)\ge 0$，所以：

$$
|z(t)|=A(t)
$$

这就是 Hilbert 变换可以表示包络的核心原理：Hilbert 变换给原实信号补出一个正交的虚部，使窄带振荡从实轴上的来回摆动变成复平面中的旋转向量；这个旋转向量的长度就是慢变化包络。

同样也可以从 Bedrosian 定理理解这一点。当 $A(t)$ 是低频信号、$\cos(2\pi f_ct)$ 是高频载波，并且二者频谱不重叠时：

$$
\mathcal{H}\{A(t)\cos(2\pi f_ct)\}
=
A(t)\mathcal{H}\{\cos(2\pi f_ct)\}
$$

由于：

$$
\mathcal{H}\{\cos(2\pi f_ct)\}=\sin(2\pi f_ct)
$$

所以：

$$
\mathcal{H}\{x(t)\}=A(t)\sin(2\pi f_ct)
$$

解析信号为：

$$
\begin{aligned}
z(t)
&=x(t)+j\mathcal{H}\{x(t)\} \\
&=A(t)\cos(2\pi f_ct)+jA(t)\sin(2\pi f_ct) \\
&=A(t)e^{j2\pi f_ct}
\end{aligned}
$$

因此仍然得到：

$$
|z(t)|=A(t)
$$

需要注意的是，这个结论依赖“单分量、窄带、包络慢变化”的条件。如果信号包含多个相互接近的分量，或者包络变化太快导致频谱和载波附近频谱严重重叠，$|z(t)|$ 仍然是解析信号的瞬时幅度，但它不一定等于直观意义上的物理包络。

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
examples/output/hilbert_transform_example.png
```
