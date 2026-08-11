# 补丁目录

- `0001-gudhi-3.13-small-n-neighbors.patch`：历史 small-N 参考补丁，基线是 GUDHI 3.13.0。它不是新 Python 库的实现基线，也不是当前 `F:\GUDHI\gudhi-devel` 开发提交上已验证的补丁。

新库的核心实现应保存在自己的源码目录中。此目录只保留历史参考或用于对照实验的补丁，不把“持续修改外部 GUDHI”作为产品架构。
