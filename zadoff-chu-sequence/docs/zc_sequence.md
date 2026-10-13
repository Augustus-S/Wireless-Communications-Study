# Zadoff-Chu 序列原理与 Python 示例解析

## 1. 概述

Zadoff-Chu sequence，简称 ZC 序列，是一类复指数序列。它最重要的特征可以概括为 CAZAC：

```text
CAZAC = Constant Amplitude Zero Auto-Correlation
```

也就是：

- Constant Amplitude：恒包络，每个采样点幅度都为 1。
- Zero Auto-Correlation：零循环自相关，除 0 移位外，与自身循环移位版本的相关结果理论上为 0。

ZC 序列还具有较好的互相关性质，并且 DFT 后幅度仍然平坦。因此它很适合用于无线通信系统中的同步、前导检测、随机接入、参考信号和信道估计等任务。

从接收机角度看，ZC 序列的价值在于“容易被相关器检测”。接收端保存一份本地 ZC 序列，然后与接收信号做相关。如果接收信号中存在对应的 ZC 序列，相关输出会出现明显峰值；峰值的位置可用于估计时延、同步点或循环移位编号。

本目录的 Python 示例会验证四件事：

1. ZC 序列样本位于单位圆上，即恒包络。
2. 循环自相关在 0 移位处有尖峰，其他移位接近 0。
3. DFT 幅度接近平坦。
4. 接收一个带噪声的循环移位 ZC 序列时，可以通过相关峰估计循环移位。

## 2. 数学定义和推导逻辑

### 2.1 序列定义

当序列长度 `N` 为奇数时，常用 ZC 序列定义为：

```text
x_u[n] = exp(-j * pi * u * n * (n + 1) / N),  n = 0, 1, ..., N-1
```

其中：

- `N` 是序列长度。
- `u` 是根索引 root index。
- `gcd(u, N) = 1`，即 `u` 和 `N` 互质。
- `j` 是虚数单位。

代码中对应的生成逻辑是：

```python
n = np.arange(length)
phase = -np.pi * root * n * (n + 1) / length
return np.exp(1j * phase)
```

当长度为偶数时，常见写法会改为：

```text
x_u[n] = exp(-j * pi * u * n^2 / N)
```

示例函数同时支持奇数长度和偶数长度：

```python
if length % 2 == 0:
    phase = -np.pi * root * n * n / length
else:
    phase = -np.pi * root * n * (n + 1) / length
```

本文档主要使用默认参数：

```python
N = 63
u = 25
```

因为 `gcd(25, 63) = 1`，满足根索引和序列长度互质的要求。

### 2.2 为什么要求 `gcd(u, N) = 1`

ZC 序列的相位项是一个二次函数：

```text
phi[n] = -pi * u * n * (n + 1) / N
```

循环相关时，本质上会比较 `x[n]` 和循环移位后的 `x[n-k]`。两个二次相位相减后，二次项抵消，剩下与 `n` 相关的一次相位项。这个一次相位项在 `n = 0 ... N-1` 上求和时，相当于一组等间隔复指数相加。

如果 `u` 与 `N` 互质，非零移位对应的复指数会完整绕单位圆分布，求和后抵消为 0。这就是零循环自相关的直观来源。

如果 `u` 与 `N` 不互质，复指数不会均匀覆盖完整周期，抵消不彻底，非零移位处会出现明显旁瓣。因此代码中直接检查：

```python
if math.gcd(root, length) != 1:
    raise ValueError("root and length must be coprime")
```

### 2.3 从相位设计理解 ZC 序列

ZC 序列不是随机序列，而是精心设计的“二次相位”序列。

普通单频复指数：

```text
s[n] = exp(j * 2*pi*f*n)
```

相邻采样点相位差是常数，因此频率固定。

ZC 序列：

```text
x[n] = exp(j * phi[n])
```

其中 `phi[n]` 是关于 `n` 的二次函数，因此相邻采样点相位差随 `n` 线性变化。这种结构类似离散 chirp。它的好处是：时域幅度恒定，同时循环相关时又能产生非常尖锐的相关峰。

## 3. 相关数学性质

### 3.1 恒包络性质

ZC 序列每个点都是复指数：

```text
x[n] = exp(j * theta[n])
```

复指数的模恒为 1：

```text
|x[n]| = 1
```

代码验证方式：

```python
magnitude = np.abs(zc)
max_magnitude_error = np.max(np.abs(magnitude - 1.0))
```

如果 `max_magnitude_error` 接近机器精度，例如 `1e-15` 量级，就说明生成结果满足恒包络。

恒包络的工程意义是：序列本身不会产生很大的幅度起伏，因此适合发射端功放处理，也不会因为局部高峰值导致明显的功率浪费。

### 3.2 循环自相关性质

普通线性相关会受到边界影响，而 ZC 序列的关键性质是循环自相关。循环移位定义为：

```text
x_k[n] = x[(n - k) mod N]
```

循环自相关为：

```text
R_xx[k] = sum_{n=0}^{N-1} x[n] * conj(x[(n - k) mod N])
```

理想 ZC 序列满足：

```text
R_xx[0] = N
R_xx[k] = 0, k != 0
```

如果做归一化：

```text
R_xx_norm[k] = R_xx[k] / N
```

则理论结果为：

```text
R_xx_norm[0] = 1
R_xx_norm[k] = 0, k != 0
```

这就是 ZC 序列适合同步和前导检测的核心原因。接收机做相关搜索时，正确位置会出现尖峰，错误位置相关值很低。

### 3.3 频域恒包络性质

ZC 序列不仅时域恒包络，DFT 后的频域幅度也接近平坦。示例中计算：

```python
spectrum = np.fft.fftshift(np.fft.fft(zc)) / np.sqrt(N)
spectrum_magnitude = np.abs(spectrum)
```

除以 `sqrt(N)` 是为了让频域幅度更便于观察。图中频谱幅度应该接近一条水平线。

这说明 ZC 序列能量在频域分布均匀，不集中在少数频点上。因此它适合在 OFDM 系统中作为同步信号、参考信号或前导序列的一部分。

### 3.4 不同根序列的互相关

不同根索引生成的 ZC 序列之间不会像自身循环自相关那样完全正交，但通常具有较低互相关。

示例中生成两个不同根序列：

```python
zc = zadoff_chu_sequence(N, root=25)
zc_other = zadoff_chu_sequence(N, root=29)
cross_corr = cyclic_correlation(zc, zc_other) / N
```

这种较低互相关有助于接收机区分多个候选序列。例如多个用户、多个小区、多个前导或多个同步候选可以使用不同根或不同循环移位构造。

### 3.5 循环移位可区分性

由于非零循环自相关为 0，同一根 ZC 序列的不同循环移位版本也可以作为不同候选序列。

接收端用本地根序列做相关：

```text
peak position -> cyclic shift
```

因此一个根序列可以通过循环移位扩展出多个可检测的前导候选。真实系统中还会考虑多径时延扩展，避免不同循环移位之间因为信道扩展而互相重叠。

## 4. 示例代码解析

示例文件：

```text
examples/zc_sequence_example.py
```

运行方式：

```bash
python3 zadoff-chu-sequence/examples/zc_sequence_example.py
```

脚本会输出数值验证结果，并保存图片：

```text
examples/output/zc_sequence_properties.png
```

### 4.1 生成 ZC 序列

核心函数如下：

```python
def zadoff_chu_sequence(length: int, root: int) -> np.ndarray:
    if length <= 0:
        raise ValueError("length must be positive")
    if math.gcd(root, length) != 1:
        raise ValueError("root and length must be coprime")

    n = np.arange(length)
    if length % 2 == 0:
        phase = -np.pi * root * n * n / length
    else:
        phase = -np.pi * root * n * (n + 1) / length
    return np.exp(1j * phase)
```

解析：

- `np.arange(length)` 生成采样索引 `n = 0 ... N-1`。
- `phase` 按 ZC 的二次相位公式生成。
- `np.exp(1j * phase)` 把相位转换成复指数序列。
- `gcd(root, length)` 检查根索引是否有效。

默认生成：

```python
length = 63
root = 25
zc = zadoff_chu_sequence(length, root)
```

### 4.2 计算循环相关

循环相关函数如下：

```python
def cyclic_correlation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if a.shape != b.shape:
        raise ValueError("sequences must have the same length")

    length = a.size
    result = np.empty(length, dtype=np.complex128)
    for shift in range(length):
        result[shift] = np.vdot(np.roll(b, shift), a)
    return result
```

其中：

- `np.roll(b, shift)` 表示对 `b` 做循环移位。
- `np.vdot()` 会对第一个输入取共轭。
- 对每个 `shift` 计算一次内积，就得到完整循环相关曲线。

当 `a == b == zc` 时：

```python
auto_corr = cyclic_correlation(zc, zc) / length
```

得到归一化循环自相关。

代码进一步计算非零移位最大旁瓣：

```python
nonzero_shifts = np.arange(1, length)
max_auto_sidelobe = float(np.max(np.abs(auto_corr[nonzero_shifts])))
```

如果 `max_auto_sidelobe` 非常接近 0，就验证了 ZC 的零循环自相关性质。

### 4.3 验证恒包络和频域平坦

恒包络验证：

```python
magnitude = np.abs(zc)
max_magnitude_error = float(np.max(np.abs(magnitude - 1.0)))
```

频域幅度：

```python
spectrum = np.fft.fftshift(np.fft.fft(zc)) / np.sqrt(length)
spectrum_magnitude = np.abs(spectrum)
```

图中第一幅子图显示复平面单位圆上的 ZC 样本点，第三幅子图显示 DFT 幅度。

### 4.4 验证不同根序列互相关

代码生成另一个根序列：

```python
other_root = 29
zc_other = zadoff_chu_sequence(length, other_root)
cross_corr = cyclic_correlation(zc, zc_other) / length
max_cross_corr = float(np.max(np.abs(cross_corr)))
```

打印 `max_cross_corr` 可以观察不同根序列之间的最大相关幅度。它不会像自相关旁瓣那样接近 0，但通常显著低于自相关主峰。

### 4.5 循环移位检测

示例模拟发送循环移位后的 ZC 序列：

```python
true_shift = 17
received = np.roll(zc, true_shift)
noise_power = 0.15
noise = noise_power * (rng.standard_normal(length) + 1j * rng.standard_normal(length))
received_noisy = received + noise
```

接收端做相关搜索：

```python
detection_corr = cyclic_correlation(received_noisy, zc) / length
estimated_shift = int(np.argmax(np.abs(detection_corr)))
```

解析：

- `true_shift = 17` 是发送端设置的循环移位。
- `received_noisy` 是加入复高斯噪声后的接收序列。
- `np.argmax(np.abs(detection_corr))` 找到最大相关峰位置。

如果 `estimated_shift == true_shift`，说明接收端正确检测到了循环移位。

### 4.6 图像输出说明

脚本输出一张包含四个子图的图片：

```text
examples/output/zc_sequence_properties.png
```

四个子图分别是：

1. ZC sequence in complex plane
   - 显示 ZC 序列点分布在单位圆上。
2. Cyclic autocorrelation
   - 显示 0 移位处有尖峰，其他移位接近 0。
3. DFT magnitude
   - 显示频域幅度接近平坦。
4. Shift detection by correlation
   - 显示接收带噪循环移位序列时，相关峰出现在真实移位位置。

## 5. 常见应用场景

### 5.1 时间同步和帧同步

接收机需要知道一帧信号从哪里开始。发送端可以在帧头放置 ZC 序列，接收端用本地 ZC 序列与接收信号做滑动相关。

当相关输出出现明显峰值时，峰值位置就是候选帧起点：

```text
received signal -> correlate with local ZC -> peak index -> frame timing
```

ZC 序列的零循环自相关特性使峰值尖锐、旁瓣低，有利于降低误同步概率。

### 5.2 随机接入前导

蜂窝系统中，终端初次接入基站时需要发送随机接入前导。ZC 序列适合做这类前导，因为：

- 多个循环移位可以构造多个前导候选。
- 接收端可以通过相关峰检测哪个前导被发送。
- 峰值位置还可以携带到达时延信息，用于上行定时调整。

LTE PRACH 就使用了基于 ZC 序列的随机接入前导设计。工程实现中会根据小区半径和多径时延扩展设置循环移位间隔，避免不同前导之间互相干扰。

### 5.3 小区搜索和同步信号

移动终端开机后需要搜索附近小区，并完成时间/频率同步。具有良好相关特性的序列适合作为同步信号。

ZC 序列或类似 CAZAC 序列可以帮助接收机：

- 判断是否存在目标同步信号。
- 估计同步信号到达时间。
- 区分不同根索引或不同小区候选。

实际标准中的同步信号可能经过截断、映射或组合，不一定直接发送完整原始 ZC 序列，但核心思想仍然是利用低旁瓣相关峰做检测。

### 5.4 参考信号和信道估计

接收机要估计无线信道，通常需要发送接收双方都已知的参考信号。ZC 序列适合做参考信号的原因是：

- 频域能量分布较均匀，有助于覆盖多个频率位置。
- 恒包络有利于发射功率稳定。
- 良好相关特性有助于从接收信号中分离参考序列。

接收端可以比较接收到的参考序列和本地参考序列，从而估计信道幅度、相位和时延结构。

### 5.5 多用户或多前导区分

同一个根 ZC 序列的不同循环移位版本可以作为不同候选；不同根索引也可以生成不同候选。这样可以构造一组相关性较低的序列集合。

应用上可以用于：

- 多个终端的随机接入前导区分。
- 多个同步候选区分。
- 多个测距或探测序列区分。

需要注意的是，真实系统不会无限制地使用任意循环移位。循环移位间隔必须大于信道最大时延扩展，否则一个用户的多径可能落到另一个循环移位候选附近，造成误检。

### 5.6 雷达、声学和测距中的相关探测

ZC 序列的低旁瓣相关特性不仅适合无线通信，也适合相关探测问题。例如：

- 发送已知 ZC 序列。
- 接收端与本地序列相关。
- 根据相关峰位置估计传播时延。

传播时延乘以传播速度即可得到距离估计。由于 ZC 序列旁瓣低，强目标旁瓣不容易掩盖弱目标或相邻路径。

### 5.7 常见注意事项

根索引必须与长度互质。如果 `gcd(u, N) != 1`，CAZAC 性质会下降。示例代码会直接抛出异常，避免生成错误序列。

ZC 的“零自相关”指的是循环自相关。使用线性相关时，边界截断会产生不同结果。

加入噪声后，非峰值位置不再精确为 0，但峰值仍应明显高于旁瓣。信噪比越高，检测越稳定。

真实通信系统中，ZC 序列可能被映射到子载波、加入循环前缀、经过信道、多径和频偏。本文示例聚焦最核心的数学性质，便于先建立直觉。

