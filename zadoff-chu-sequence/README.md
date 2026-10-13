# Zadoff-Chu Sequence Study

本目录用于学习 Zadoff-Chu sequence（ZC 序列）的数学原理、相关特性和 Python 仿真验证。

ZC 序列是一类恒包络、零循环自相关（CAZAC, Constant Amplitude Zero Auto-Correlation）复序列。它常用于无线通信中的同步、随机接入前导、参考信号和序列检测等场景。典型系统包括 LTE/NR 中的同步和随机接入相关设计。

## 目录结构

```text
zadoff-chu-sequence/
  README.md
  requirements.txt
  docs/
    zc_sequence.md
  examples/
    zc_sequence_example.py
    output/
```

## 快速运行

安装依赖：

```bash
pip install -r zadoff-chu-sequence/requirements.txt
```

运行示例：

```bash
python3 zadoff-chu-sequence/examples/zc_sequence_example.py
```

脚本会：

- 生成长度为 `N=63`、根索引 `u=25` 的 ZC 序列。
- 验证恒包络特性。
- 计算循环自相关，展示非零循环移位处接近 0。
- 计算不同根序列之间的循环互相关。
- 模拟接收端通过相关峰检测循环移位。
- 将图保存到 `zadoff-chu-sequence/examples/output/zc_sequence_properties.png`。

## 文档

- [docs/zc_sequence.md](docs/zc_sequence.md): 结合 Python 代码逐步讲解 ZC 序列定义、CAZAC 特性、频域性质、互相关和循环移位检测。

