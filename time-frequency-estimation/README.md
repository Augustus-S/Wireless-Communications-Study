# Time-Frequency Estimation Examples

本项目给出 4 个 Python 时频估计示例，并配套说明数学原理：

- `examples/stft_example.py`: 短时傅里叶变换 STFT
- `examples/wavelet_transform_example.py`: Morlet 连续小波变换 CWT
- `examples/hilbert_transform_example.py`: Hilbert 变换与解析信号
- `examples/instantaneous_frequency_example.py`: 基于解析相位的瞬时频率估计

## 快速运行

```bash
python3 examples/stft_example.py
python3 examples/wavelet_transform_example.py
python3 examples/hilbert_transform_example.py
python3 examples/instantaneous_frequency_example.py
```

脚本会将图片保存到 `outputs/`。

## 文档

- `docs/stft.md`
- `docs/wavelet_transform.md`
- `docs/hilbert_transform.md`
- `docs/instantaneous_frequency.md`
