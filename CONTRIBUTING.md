# 为 Topp 贡献

[English](CONTRIBUTING.en.md) · [开发指南](docs/DEVELOPMENT.md)

欢迎提交正确性、性能、跨平台构建、测试、示例和文档改进。

1. 保持修改范围小，并说明它属于公开 Python API、稳定默认内核还是实验策略。
2. 公共行为变更必须同时更新中英文 API/Usage、类型提示和 Python 测试。
3. 内核优化必须先通过 C++ reference、Python 契约和适用的外部 oracle 差分，再报告随机执行顺序的 median/p95。
4. 不提交 `build/`、`.obj`、wheel、原始 benchmark 结果、虚拟环境或第三方源码树。
5. 性能结论必须注明提交、编译器、CPU、数据分布、规模、重复次数与正确性边界。

提交前至少运行：

```powershell
py -m pytest tests/python
cmd.exe /d /c scripts\build-kernel.cmd
build\manual\bottleneck_core_tests.exe
build\manual\wasserstein_core_tests.exe
```
