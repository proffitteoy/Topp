# 内核研究归档

本目录保存已经完成的内核实验结论和第三方对照补丁，不是普通用户使用文档。Topp 1.0 的生产路径、保留维护入口和已知边界先看[内核最终状态](FINAL_STATE.md)。

## 阅读边界

- `evidence/` 记录特定提交、机器、编译器和输入上的实验；数字只在文档注明的边界内成立。
- `patches/` 保存用于对照的第三方历史补丁；不进入 Topp 运行时，也不是当前实现基线。
- 已完成或被否定的长篇方案草稿不再保留在活动树中；需要追溯时使用 Git 历史。

发布级用法与契约以根目录 README、`docs/USAGE.md`、`docs/API.md` 和测试为准。研究文档是按日期冻结的证据；若与当前发布文档冲突，以当前发布文档和测试为准。

## 冻结结论

- [Topp 1.0 内核最终状态](FINAL_STATE.md)
- [Wasserstein 内核实验与主线冻结](evidence/WASSERSTEIN_KERNEL.md)
- [Bottleneck adaptive router 统一消融](evidence/PHASE2_ROUTER_REPORT.md)
- [Bottleneck 第一阶段 exact 内核](evidence/PHASE1_KERNEL_REPORT.md)

## 分项证据

- [初始化基线](evidence/BASELINE.md)
- [Bottleneck 实验登记](evidence/KERNEL_EXPERIMENTS.md)
- [multiplicity capacity](evidence/PHASE2_MULTIPLICITY_REPORT.md)
- [mandatory partial flow](evidence/PHASE2_MANDATORY_FLOW_REPORT.md)
- [geometry 与 component](evidence/PHASE2_GEOMETRY_COMPONENT_REPORT.md)
- [exact threshold search](evidence/PHASE2_SEARCH_REPORT.md)
- [one-to-many preparation 与 batch reuse](evidence/PHASE2_BATCH_REPORT.md)
- [Wasserstein assist 负结果](evidence/PHASE2_WASSERSTEIN_ASSIST_REPORT.md)
- [历史补丁说明](patches/README.md)
