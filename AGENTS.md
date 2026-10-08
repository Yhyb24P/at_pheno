# 项目开发约定

- 阅读 README.md、preprocessing/README.md 与 docs/研究设计.md；使用 preprocessing/processed 中的数据。
- Python 环境使用 Miniconda、environment.yml 和本地 .venv。
- 优先修改现有实现，减少独立脚本、重复文档和目录层级。
- 数据处理文件集中在 preprocessing；已完成的检查不重复执行，除非输入或相关实现发生变化。
- 保留样本对齐、缺失值语义及训练/测试隔离；频率筛选和填补在训练分区内完成。
- 重要数值改动运行现有针对性测试；文档和文件整理不新增测试。通过后继续工作。
- 原始数据和可用产物保持稳定，不入 Git；运行输出写入独立目录。
