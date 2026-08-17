# 研究归档

本目录保存内核研究证据、历史方案和实验补丁，不是普通用户使用文档。

## 状态

- `evidence/`：在特定提交、机器、编译器和输入上获得的实验记录；结论只在文档注明的边界内成立。
- `proposals/`：历史或未来设计方案；内容可能尚未实现，也可能已被后续实验否定。
- `patches/`：用于对照的历史第三方补丁；不进入 Topp 运行时，也不是当前实现基线。

发布级用法和契约以根目录 README、`docs/USAGE*` 与 `docs/API*` 为准。若研究文档与发布级文档冲突，以发布级文档和测试为准。

## Evidence

- [初始化基线](evidence/BASELINE.md)
- [Bottleneck 实验登记](evidence/KERNEL_EXPERIMENTS.md)
- [Bottleneck 第一阶段报告](evidence/PHASE1_KERNEL_REPORT.md)
- [Bottleneck 第二阶段 multiplicity 报告](evidence/PHASE2_MULTIPLICITY_REPORT.md)
- [Bottleneck 第二阶段 mandatory partial flow 报告](evidence/PHASE2_MANDATORY_FLOW_REPORT.md)
- [Bottleneck 第二阶段 geometry/component 报告](evidence/PHASE2_GEOMETRY_COMPONENT_REPORT.md)
- [Bottleneck 第二阶段 search strategy 报告](evidence/PHASE2_SEARCH_REPORT.md)
- [Bottleneck 第二阶段 batch reuse 报告](evidence/PHASE2_BATCH_REPORT.md)
- [Bottleneck 第二阶段 adaptive router 报告](evidence/PHASE2_ROUTER_REPORT.md)
- [Bottleneck 第二阶段 Wasserstein assist 报告](evidence/PHASE2_WASSERSTEIN_ASSIST_REPORT.md)
- [Wasserstein 内核实验记录](evidence/WASSERSTEIN_KERNEL.md)

## Proposals

- [项目范围历史草案](proposals/PROJECT_SCOPE.md)
- [Large-N 特化](proposals/large-N特化.md)
- [Wasserstein 优化方案](proposals/wasserstein优化方案.md)
- [Wasserstein 第二阶段](proposals/wasserstein第二阶段.md)
- [Bottleneck 优化方向](proposals/优化方向.md)
- [Bottleneck 第二阶段](proposals/bottleneck第二阶段.md)
- [历史发布流程建议](proposals/发布流程.md)
- [历史 Python API 长篇草案](proposals/接口文档.md)
