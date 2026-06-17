# 瞬时频率估计

## 1. 方法目标

瞬时频率估计试图直接估计信号在每一个时刻的局部频率。对于单分量窄带信号，可以先通过 Hilbert 变换构造解析信号，再对解析相位求导。

这种方法适合单分量、信噪比较高、相位变化连续的信号；对强噪声、多分量叠加和幅度过零比较敏感。

## 2. 数学原理

解析信号写作：

$$
z(t)=x(t)+j\mathcal{H}\{x(t)\}=A(t)e^{j\phi(t)}
$$

其中 $\phi(t)$ 是瞬时相位。瞬时角频率定义为：

$$
\omega_i(t)=\frac{d\phi(t)}{dt}
$$

瞬时频率为：

$$
f_i(t)=\frac{1}{2\pi}\frac{d\phi(t)}{dt}
$$

对于线性调频信号：

$$
x(t)=\cos\left(2\pi\left(f_0t+\frac{1}{2}kt^2\right)\right)
$$

其相位为：

$$
\phi(t)=2\pi\left(f_0t+\frac{1}{2}kt^2\right)
$$

因此瞬时频率是：

$$
f_i(t)=\frac{1}{2\pi}\frac{d\phi(t)}{dt}=f_0+kt
$$

离散时间中，如果采样率为 $f_s$，相邻采样点间隔为：

$$
\Delta t=\frac{1}{f_s}
$$

则可以用差分近似：

$$
\hat{f}_i[n]\approx\frac{f_s}{2\pi}\left(\phi[n]-\phi[n-1]\right)
$$

代码中使用 `np.gradient`，相当于中心差分形式：

$$
\hat{f}_i[n]\approx\frac{f_s}{2\pi}\frac{\phi[n+1]-\phi[n-1]}{2}
$$

## 3. 与代码的对应关系

示例脚本：`examples/instantaneous_frequency_example.py`

`make_noisy_chirp(fs, duration, rng)` 构造带噪线性调频信号：

$$
x(t)=\cos\left(2\pi\left(f_0t+\frac{1}{2}kt^2\right)\right)+\eta(t)
$$

其中 $\eta(t)$ 是随机噪声。函数同时返回真实瞬时频率：

$$
f_i(t)=f_0+kt
$$

`estimate_instantaneous_frequency(x, fs, smoothing_length)` 首先构造解析信号：

```python
z = signal.hilbert(x)
```

对应：

$$
z(t)=x(t)+j\mathcal{H}\{x(t)\}
$$

然后计算并展开相位：

```python
phase = np.unwrap(np.angle(z))
```

对应：

$$
\phi(t)=\mathrm{unwrap}(\arg z(t))
$$

接着对相位求导：

```python
instantaneous_frequency = np.gradient(phase) * fs / (2.0 * np.pi)
```

对应：

$$
\hat{f}_i(t)=\frac{1}{2\pi}\frac{d\phi(t)}{dt}
$$

由于噪声会被相位求导放大，代码使用中值滤波：

```python
smoothed = signal.medfilt(instantaneous_frequency, kernel_size=smoothing_length)
```

它不改变瞬时频率定义，只是对估计轨迹做后处理，使曲线更稳定。

## 4. 与 STFT 谱脊估计的区别

STFT 谱脊估计通常从二维能量图中取最大值：

$$
\hat{f}(t)=\arg\max_f P(t,f)
$$

Hilbert 瞬时频率估计直接从相位求导：

$$
\hat{f}_i(t)=\frac{1}{2\pi}\frac{d\phi(t)}{dt}
$$

前者更像从时频图中提取轨迹；后者更像从信号本身的解析相位中提取局部频率。

## 5. 运行方式

```bash
python3 examples/instantaneous_frequency_example.py
```

输出图片：

```text
outputs/instantaneous_frequency_example.png
```
