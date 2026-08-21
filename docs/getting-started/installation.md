# 安装

[English](https://proffitteoy.github.io/Topp/en/getting-started/installation.html)

## Windows 与 Linux wheel

Topp 为 Windows x64、Linux x86_64 和 CPython 3.10–3.14 发布 wheel：

```console
py -m pip install topp
```

确认安装版本：

```pycon
>>> import topp
>>> topp.__version__
'1.0.1'
```

NumPy 1.23 或更高版本会自动安装，也是 Topp 唯一的运行时依赖。

## 虚拟环境

建议使用隔离环境：

```console
py -m venv .venv
.venv\Scripts\activate
py -m pip install topp
```

## 从源码构建

源码构建需要 CMake 3.24 或更高版本，以及支持 C++20 的编译器：

```console
git clone https://github.com/proffitteoy/Topp.git
cd Topp
py -m pip install -v .
```

Linux 已纳入发布 CI；macOS 尚未纳入。在 macOS 成功完成本地编译，并不代表它属于当前公开支持范围。

## 冒烟测试

```python
import topp

assert topp.bottleneck_distance([[0.0, 1.0]], []) == 0.5
print("Topp 已就绪")
```

预期输出：

```text
Topp 已就绪
```
