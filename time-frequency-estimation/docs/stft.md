# STFT 短时傅里叶变换

## 1. 方法目标

短时傅里叶变换用于估计信号的频率成分如何随时间变化。它把长信号切成一系列局部时间窗口，并对每个窗口分别做傅里叶变换，因此得到二维时频表示。

对于非平稳信号，普通傅里叶变换只能给出全局频谱：

$$
X(f)=\int_{-\infty}^{+\infty}x(t)e^{-j2\pi ft}\,dt
$$

它不能直接回答“某个频率在什么时候出现”。STFT 的输出则是：

$$
X(\tau,f)
$$

其中 $\tau$ 表示窗口中心时间，$f$ 表示频率。

## 2. 数学原理

连续时间 STFT 定义为：

$$
X(\tau,f)=\int_{-\infty}^{+\infty}x(t)w(t-\tau)e^{-j2\pi ft}\,dt
$$

其中：

- $x(t)$ 是待分析信号。
- $w(t-\tau)$ 是以 $\tau$ 为中心的窗函数。
- $e^{-j2\pi ft}$ 是傅里叶基函数。
- $X(\tau,f)$ 是时间 $\tau$、频率 $f$ 处的复数谱系数。

常见时频谱图使用功率谱：

$$
P(\tau,f)=|X(\tau,f)|^2
$$

代码中使用分贝显示：

$$
P_{\mathrm{dB}}(\tau,f)=10\log_{10}(P(\tau,f)+\epsilon)
$$

其中 $\epsilon$ 是很小的正数，用于避免 $\log 0$。

STFT 的核心折中是时间分辨率和频率分辨率：

$$
\Delta t \cdot \Delta f \gtrsim C
$$

窗口越短，时间定位越准，但频率分辨率越差；窗口越长，频率分辨率越高，但时间定位越模糊。

## 3. 与代码的对应关系

示例脚本：`examples/stft_example.py`

`make_chirp(fs, duration)` 构造线性调频信号：

$$
x(t)=\cos\left(2\pi\left(f_0t+\frac{1}{2}kt^2\right)\right)
$$

其瞬时频率为：

$$
f_i(t)=f_0+kt
$$

这里的真实瞬时频率用于和 STFT 谱脊估计结果对比。

`compute_stft(x, fs, window_length, overlap)` 调用：

```python
signal.stft(
    x,
    fs=fs,
    window="hann",
    nperseg=window_length,
    noverlap=overlap,
    nfft=512,
    boundary=None,
)
```

对应公式中的窗函数 $w(t-\tau)$ 和局部傅里叶变换。`nperseg` 控制窗口长度，`noverlap` 控制相邻窗口重叠，`nfft` 控制频率采样点数。

`power = np.abs(spectrum) ** 2` 对应：

$$
P(\tau,f)=|X(\tau,f)|^2
$$

`estimate_ridge(frequencies, power)` 使用每个时间帧中功率最大的频率作为频率轨迹估计：

$$
\hat{f}(\tau)=\arg\max_f P(\tau,f)
$$

这是一个最简单的时频估计方法，适合单分量或主分量明显的信号。

## 4. 运行方式

```bash
python3 examples/stft_example.py
```

输出图片：

```text
examples/output/stft_example.png
```
