## 编译

### 仅编译 ddsim 库（不含可执行文件）

```sh
cmake -S . -B /tmp/ddsim -DBUILD_DDSIM=ON -DGIT_SUBMODULE=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/ddsim --target ddsim -j$(nproc)
```

### 编译 ddsim 库 + 可执行文件（simple, noise_aware, primebases, vectors, benchmark）

```sh
cmake -S . -B /tmp/ddsim -DBUILD_DDSIM=ON -DBUILD_DDSIM_APPS=ON -DGIT_SUBMODULE=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/ddsim -j$(nproc)
```

### 编译全部项目（qfr + ddsim + 测试）

```sh
cmake -S . -B /tmp/ddsim -DBUILD_DDSIM=ON -DBUILD_DDSIM_APPS=ON -DBUILD_QFR_TESTS=ON -DGIT_SUBMODULE=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/ddsim -j$(nproc)
```

### 生成的 ddsim 可执行文件

| 文件 | 说明 |
|------|------|
| `ddsim_simple` | 基本模拟器前端 |
| `ddsim_noise_aware` | 噪声感知模拟器 |
| `ddsim_primebases` | Shor 算法素因数分解 |
| `ddsim_vectors` | 向量输出工具 |
| `ddsim_benchmark` | 性能基准测试 |
| `ddsim_test` | 单元测试 |