# 清单与哈希

- `TREE.txt`：压缩前的文件树。
- `FILE_MANIFEST_SHA256.tsv`：相对路径、文件字节数和 SHA256。
- `PACKAGE_VALIDATION.json`：目录、release、图和 docking 的交接前验证结果。

`FILE_MANIFEST_SHA256.tsv` 不包含它自身，以避免循环哈希。压缩包本身的 SHA256 放在压缩包旁边的 `finale.zip.sha256.txt`。
