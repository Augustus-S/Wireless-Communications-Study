# 小波变换

## 1. 方法目标

小波变换用于多尺度时频分析。与 STFT 使用固定窗口不同，小波变换通过尺度参数改变分析窗口：

- 高频信号使用短窗口，时间定位更好。
- 低频信号使用长窗口，频率定位更好。

因此小波变换适合分析瞬态、脉冲、突发和多尺度信号。

## 2. 数学原理

连续小波变换定义为：

$$
W_x(a,b)=\frac{1}{\sqrt{|a|}}\int_{-\infty}^{+\infty}x(t)\psi^*\left(\frac{t-b}{a}\right)\,dt
$$

其中：

- $\psi(t)$ 是母小波。
- $a$ 是尺度参数。
- $b$ 是时间平移参数。
- $^*$ 表示复共轭。
- $W_x(a,b)$ 是尺度 $a$、时间 $b$ 处的小波系数。

尺度 $a$ 与频率成反比。尺度越小，小波越窄，对应更高频；尺度越大，小波越宽，对应更低频。

本示例使用复 Morlet 小波：

$$
\psi(t)=\pi^{-1/4}e^{j\omega_0 t}e^{-t^2/2}
$$

尺度化后：

$$
\psi_a(t)=\frac{1}{\sqrt{a}}\psi\left(\frac{t}{a}\right)
$$

对于 Morlet 小波，尺度到伪频率的近似关系为：

$$
f_a \approx \frac{\omega_0 f_s}{2\pi a}
$$

其中 $f_s$ 是采样率，$\omega_0$ 是 Morlet 小波的中心角频率参数。

小波时频能量图通常取：

$$
P(a,b)=|W_x(a,b)|^2
$$

如果把尺度换算成伪频率，也可写为：

$$
P(f,b)=|W_x(f,b)|^2
$$

## 3. 与代码的对应关系

示例脚本：`examples/wavelet_transform_example.py`

`make_multiscale_signal(fs, duration)` 构造包含低频成分、高频短时突发和高斯脉冲的信号，用于展示小波对多尺度结构的响应。

`morlet_cwt(x, fs, frequencies, omega0)` 实现连续小波变换。代码先根据目标伪频率计算尺度：

```python
scales = omega0 * fs / (2.0 * np.pi * frequencies)
```

对应公式：

$$
a=\frac{\omega_0 f_s}{2\pi f_a}
$$

每个尺度上生成离散 Morlet 小波：

```python
wavelet = np.pi ** (-0.25) * np.exp(1j * omega0 * n / scale)
wavelet *= np.exp(-0.5 * (n / scale) ** 2) / np.sqrt(scale)
```

对应：

$$
\psi_a[n]=\frac{1}{\sqrt{a}}\pi^{-1/4}
e^{j\omega_0 n/a}e^{-\frac{1}{2}(n/a)^2}
$$

然后用卷积计算小波系数：

```python
coefficients[row] = signal.fftconvolve(x, np.conj(wavelet[::-1]), mode="same")
```

这对应连续公式中的内积：

$$
W_x(a,b)=\langle x,\psi_{a,b}\rangle
$$

最后：

```python
power = np.abs(coefficients) ** 2
```

对应小波能量：

$$
P(a,b)=|W_x(a,b)|^2
$$

## 4. 运行方式

```bash
python3 examples/wavelet_transform_example.py
```

输出图片：

```text
outputs/wavelet_transform_example.png
```
