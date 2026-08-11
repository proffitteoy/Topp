# 初始化基线

记录日期：2026-08-10

这份文件只记录初始化时已经核实的事实。项目目标是独立 Python 库；外部历史性能数据只用于 oracle、参考实现和基准语料，不代表当前仓库已经复现，也不构成产品依赖。

## 当前仓库

- 路径：`F:\bottleneck`
- Git：仓库已存在，但初始化检查时没有提交，分支名为 `master`
- 初始化前除 `.git/` 外没有项目文件
- 本次没有安装依赖、编译 GUDHI 或修改外部源码树

## 外部 GUDHI 源码树

- 路径：`F:\GUDHI\gudhi-devel`
- 当前分支：`agent/reader-utilities-invalid-paths`
- 当前提交：`d645db79bb8bb38eb2c3020205f33a45a5f2c18a`
- `git describe`：`tags/gudhi-release-3.13.0-5-gd645db79b`
- `pyproject.toml` 声明版本：`3.14.0a0`
- 远端：个人 fork 为 `origin`，GUDHI 官方仓库为 `upstream`
- 工作树有未跟踪目录 `.vs/` 和 `out/`；后续操作必须保留它们
- Hera 头文件 `ext/hera/include/hera/bottleneck.h` 存在

当前提交并不是干净的 GUDHI 3.13.0 标签。因此，基于 3.13.0 制作的历史补丁不能被假定为可直接用于当前开发提交；选定目标版本后需要重新检查调用路径、上下文和测试。

初始化时对当前提交运行了只读的 `git apply --check`，历史补丁的文本上下文检查通过。这只说明补丁目前能够套用，不代表它已经在该提交上编译、通过测试或保持性能收益。

## 已发现的构建环境

外部源码树已有 `out/build/x64-Debug` 配置缓存：

- 生成器：Ninja
- 编译器：Visual Studio 2022 Community / MSVC 14.44
- 构建类型：Debug
- `WITH_GUDHI_TEST=ON`
- `WITH_GUDHI_BENCHMARK=OFF`
- 缓存中的 `CGAL_DIR` 和 `TBB_DIR` 为 `NOTFOUND`
- Boost 指向 `F:/GUDHI/vcpkg/installed/x64-windows/share/boost`

系统可发现 CMake 3.29.2 和 Ninja，但它们来自 MiKTeX 附带的工具目录；初始化时没有把该组合认定为正式 Release 基准工具链。后续源码编译前应创建独立的 Release 构建目录并记录完整依赖版本与编译选项。

## 历史 small-N 参考实验

历史工作在 GUDHI 3.13.0 上修改 `Neighbors_finder`：当 `Persistence_graph::size() <= 256` 时使用连续数组线性扫描和 swap-pop 删除，更大的图继续使用 CGAL kd-tree。

历史记录报告：

| 数据组 | 官方 3.13.0 | small-N | 内核/Python 入口加速 |
|---|---:|---:|---:|
| 较大图组，128 对，点数中位数 120 | 559.511 μs | 366.345 μs | 1.53x |
| 较小图组，128 对，点数中位数 18 | 144.535 μs | 102.879 μs | 1.40x |

当时还记录了两组距离二进制摘要一致、250 组随机差分无差异以及 8 项 GUDHI Bottleneck 测试通过。这些结果没有在当前仓库复现。

历史 wheel 为 CPython 3.12 / Windows x64 专用，不能作为新 Python 库的发行产物，也不能作为其他 ABI、平台或当前开发分支的验证结果。补丁副本保存在 `patches/0001-gudhi-3.13-small-n-neighbors.patch`。

## 后续方案需要明确的事项

- 用哪个固定 GUDHI 版本作为 `e=0` oracle 与性能基线
- 新库首版只提供 exact，还是也公开 approximate 接口
- 目标持久图的点数、有限/无限点、重复点和数值范围分布
- exact 语义要求是相同数学值还是浮点逐位一致
- 成功指标覆盖 C++ 内核、Python 单次调用和批量计算；具体下游应用另行报告
- 包名、平台/Python 版本矩阵与构建后端
