# MOM 库耦合到 Ansys Fluent UDF 的可行性评估

> 评估对象: acuoci/mom (C++20/23 矩方法库)
> 目标: 通过 UDF (User-Defined Functions) 将 MOM 颗粒源项接入 Fluent CFD 求解器
> 依据: Ansys Fluent UDF Manual (官方 v12.0, 架构沿用至今) + mom 仓库架构

---

## 结论 (TL;DR)

**可行, 但需要显著的集成工程, 不是即插即用。**

mom 是纯 C++20/23 模板库, Fluent UDF 构建链默认是 C89 模式 (`cl /c /Za`)。要打通, 必须在二者之间架一层 **C++ 桥接层**, 并改写 UDF makefile 以:
1. 用 C++20 编译桥接 `.cpp`
2. 链接预编译的 MOM 静态库 (`MOM.lib` / `libMOM.a`)
3. 用 `extern "C"` 暴露 UDF 入口 (Fluent 加载器只认 C 符号)
4. 处理并行分区下的每分区模型实例生命周期

不存在原则性障碍 — Fluent 官方明确支持链接预编译目标文件/外部库 (UDF Manual §node135), C++ 也已多次被社区验证可行。

---

## 1. 官方 UDF 机制关键事实

### 1.1 DEFINE_SOURCE — 单元源项钩子

```c
DEFINE_SOURCE(name, c, t, dS, eqn)
// 返回 real, 对单个单元 c 计算一个源项值
// dS[eqn] = dS/dphi, 可选, 用于隐式线性化 S = A + B·phi
// 单位 = 生成率/体积 (如 kg/m³·s)
// 仅编译型 UDF, 非解释型
```

- 求解器自动按单元线程循环调用, UDF 本身**不必**循环单元
- 适用方程: 质量、动量、k、ε、能量(含固体区)、组分质量分数、P1 辐射、**UDS 输运**、颗粒温度
- **不适用于**: 离散坐标 (DO) 辐射

mom 返回的 `std::span<const double>` 源项 (NEq 维) 可直接映射到: (a) UDS 输运方程源项 (每个矩一个 UDS), 或 (b) 组分质量分数源项 (碳烟/气相物种)。

### 1.2 UDF 库构建 (Windows)

- 编辑每个架构/版本文件夹下的 `user_nt.udf`:
  - `SOURCES = $(SRC)myudf.c ...`
  - `VERSION = 3ddp` 等
  - `PARALLEL_NODE = none|smpi|vmpi|net`
- 在 VS 命令提示符下运行 `nmake`
- 编译器标志 (实测日志): `cl /c /Za /DUDF_EXPORTING -I<fluent>/src ...`
- **`/Za` = 禁用语言扩展 (严格 ANSI C)** ← 这是 C 模式, 与 C++20 冲突的关键
- 产物: `libudf.dll`, 链接 `fl1209s.lib`

### 1.3 链接预编译目标文件 (官方支持外部库) — 关键

UDF Manual §node135 明确给出链接外部预编译 `.o`/`.obj` 的流程 (以 FORTRAN 为例):

- **Windows**: 拷贝 `.obj` 到各架构/版本目录, 编辑 `user_nt.udf`, 设置 `USER_OBJECTS = myobj1.obj myobj2.obj`
- **Linux**: 拷贝 `.o` 到各架构目录, 编辑 makefile `USER_OBJECTS = myobject1.o myobject2.o` (空格分隔)
- 官方要求: 目标文件需用与 Fluent **相似的标志**编译

这是链接 MOM 静态库的官方认可机制。

### 1.4 链接静态库 (社区验证)

- cfd-online 论坛: 需修改 makefile 在链接阶段加入 `project1.lib`
- ResearchGate: 在 makefile 的链接命令中包含 lib 名称
- 必须 `extern "C"` 包裹 UDF 入口 (`DEFINE_*` 宏展开为 C 链接符号, Fluent 加载器按 C 符号查找); **声明和定义都要 extern C**

---

## 2. mom 库相关特性回顾

| 特性 | 值 | 对 UDF 集成的影响 |
|------|----|------------------|
| 语言标准 | C++20 (核心), C++23 (字典模式) | 必须改 makefile 用 `/std:c++20` |
| 模板/CRTP | 重度, 零虚函数 | 静态库已显式实例化 `BasicThermoData`, 链接即可 |
| `extern template` | 启用 (MOM_COMPILED_LIBRARY) | 好: 编译时间可控 |
| Eigen 3.4 | header-only | 无 ABI 风险, 只需 include 路径 |
| 入口 API | `ComputeCell(model,T,P,Y,mu,moments) noexcept` | 干净的单单元调用, 与 DEFINE_SOURCE 同构 |
| 返回类型 | `std::span<const double>` | C++ 类型, 必须在桥接层内部消费, 转成 `real` |
| `AnyMomentMethod` | `std::variant` + `std::visit` | 每分区一个实例, 状态可变 |
| 依赖 | 仅 stdlib + Eigen | 轻量, 无运行时依赖 |
| 构建 | CMake, 产物 `MOM.lib`/`libMOM.a` | 标准 MSVC/GCC 静态库 |

---

## 3. 集成挑战与对策

### 挑战 1: C++20 vs UDF 默认 C 编译模式

**问题**: Fluent `user_nt.udf` 用 `cl /c /Za` (C 模式, 严格 ANSI)。MOM 需要 C++20, 需 `/std:c++20 /TP` (或源文件用 `.cpp` 扩展名让 cl 自动选 C++)。

**对策**:
- UDF 入口文件用 `.cpp` 扩展名 (cl 按扩展名选语言)
- 改写 `user_nt.udf`, 在编译器选项中加入 `/std:c++20 /EHsc`
- 去掉 `/Za` (C++ 不接受 `/Za` 的全部语义; 改用 `/permissive-` 严格标准)
- `udf.h` 本身有 `#ifdef __cplusplus` 守护, 可在 C++ 模式下包含

### 挑战 2: 名称修饰 (Name Mangling)

**问题**: `DEFINE_SOURCE` 等宏展开为 Fluent 加载器按 **C 符号**查找的入口。C++ 默认名称修饰会导致符号找不到。

**对策**:
- 桥接 `.cpp` 中所有 `DEFINE_*` 宏外层用 `extern "C"` 包裹
- `udf.h` 内部已对宏做了 `extern "C"` 守护 (若没有, 在桥接文件顶部 `extern "C" { #include "udf.h" }`)
- 对外暴露的辅助函数 (供其他 UDF 调用) 也用 `extern "C"`

### 挑战 3: 模型实例生命周期 (无自然初始化点)

**问题**: `DEFINE_SOURCE` 按单元、按方程被调用, 没有自然的"创建模型"时机。`AnyMomentMethod` 是带状态的 `std::variant`, 不能每次调用都重建。

**对策** (标准 UDF 模式):
- `DEFINE_INIT` 或 `DEFINE_ADJUST` 中创建/更新模型实例
- 句柄存入全局指针 (静态/全局变量) 或 Fluent 用户定义内存 (`C_UDMI`)
  - 全局静态指针: 简单, 但并行下需每节点一份
  - `C_UDMI`: 每单元存储, 过度; 实际用全局静态 + 分区隔离更合理
- `DEFINE_SOURCE` 内取句柄 → `SetState` → `SetMoments` → `ComputeCell` → 返回源项

### 挑战 4: 并行分区与线程安全

**问题**: Fluent 并行把网格分到多个计算节点; 每个节点是独立进程。MOM 模型对象在 `Compute` 期间持有可变状态 (`moments_`, `omega_gas_`), 不能跨线程共享。

**对策**:
- 每个 compute node 持有独立模型实例 (静态全局即可, 因为 Fluent 计算节点是独立进程, 不共享地址空间)
- 若同一节点内多线程 (Fluent 共享内存并行), 需 `thread_local` 或每线程实例
- Fluent 提供 `node_host_to_node_*` 通信; 通常模型配置在 host 广播到所有 nodes

### 挑战 5: ABI / MSVC 版本对齐

**问题**: Fluent 用其捆绑的 MSVC 工具链构建 `libudf`。MOM 静态库必须用**相同**编译器版本 + 相同 CRT 配置 (`/MD` vs `/MT`, `_DEBUG` vs `_RELEASE`) 编译, 否则链接/运行时崩溃。

**对策**:
- 查 Fluent 安装目录的编译器版本 (通常捆绑 VS 2019/2022)
- 用同版本 MSVC 构建 MOM (CMake + VS 生成器, `-G "Visual Studio 17 2022" -A x64`)
- CRT 链接方式对齐: Fluent 默认 `/MD` (多线程 DLL), MOM 构建时设 `MSVC_RUNTIME_LIBRARY=MultiThreadedDLL`
- Eigen header-only 无 ABI 问题

### 挑战 6: 源项到 UDS 方程的映射

**问题**: mom 返回 `NEq` 维源项 (HMOM=4, ThreeEquations=3, 等)。Fluent 中每个矩需对应一个 UDS 标量方程, `DEFINE_SOURCE` 按方程逐个注册。

**对策**:
- 为每个矩注册一个 UDS (Fluent 最多支持 N 个 UDS, N 通常 ≥ 4)
- 每个矩写一个 `DEFINE_SOURCE` (或一个 UDF 内用 `eqn` 区分)
- 一次 `ComputeCell` 结果缓存 (在 `DEFINE_ADJUST` 算完存 `C_UDMI`), 各 `DEFINE_SOURCE` 直接读 — 避免重复计算
- 单位换算: mom 矩单位 (mol/m³ 等) → Fluent UDS 标量单位, 注意一致性

---

## 4. 推荐集成架构

```
┌─────────────────────────────────────────────────────┐
│  libudf.dll  (Fluent 加载)                           │
│                                                      │
│  ┌──────────────────────────────────────────────┐   │
│  │  udf_bridge.cpp  (C++20, extern "C" 入口)     │   │
│  │                                                │   │
│  │  DEFINE_INIT    → 创建模型实例 (每分区)         │   │
│  │  DEFINE_ADJUST  → 每步: 算源项存 C_UDMI        │   │
│  │  DEFINE_SOURCE  → 读 C_UDMI 返回 (每矩一个)    │   │
│  │  DEFINE_RW      → 读写 restart (模型状态)      │   │
│  └──────────────┬───────────────────────────────┘   │
│                  │ 调用 C++ API                       │
│  ┌──────────────▼───────────────────────────────┐   │
│  │  MOM.lib  (预编译静态库, C++20, MSVC)          │   │
│  │  + Eigen 3.4 (header-only)                    │   │
│  │  ComputeCell / ForEachCell / GetSources        │   │
│  └──────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────┘
```

---

## 5. 七步集成计划

### 步骤 1: 构建 MOM 静态库 (匹配 Fluent MSVC)

```bat
cd E:\HMOM\mom
cmake -B build -G "Visual Studio 17 2022" -A x64 ^
  -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL ^
  -DMOM_BUILD_TESTS=OFF -DMOM_BUILD_SHARED=OFF
cmake --build build --config Release
:: 产物: build\MOM\Release\MOM.lib, 以及 Eigen 头路径
```

记录: `MOM.lib` 路径、Eigen include 路径、MSVC 版本。

### 步骤 2: 确定 Fluent UDF 构建环境

- 找到 Fluent 安装的 `src/`、`lib/` 目录
- 确认捆绑编译器版本 (运行 Fluent 自带的 UDF 编译一次, 看日志的 `cl` 版本)
- 选定目标架构文件夹 (如 `3ddp_node`)

### 步骤 3: 编写 C++ 桥接 UDF

`udf_bridge.cpp` 骨架:

```cpp
// C++20, 包含 Fluent udf.h + MOM 头
extern "C" {
#include "udf.h"      // Fluent UDF 宏, C 链接
}

#include "MOM/MOM.hpp"
#include <memory>
#include <new>

namespace {
// 每分区一个模型实例 (Fluent compute node = 独立进程)
struct ModelState {
    MOM::BasicThermoData thermo;
    MOM::AnyMomentMethod<MOM::BasicThermoData> model;
    // ... 配置参数
};

std::unique_ptr<ModelState> g_state;  // 进程级单例

constexpr int kUdmSootMomentSrcBase = 0;  // C_UDMI 起始槽
constexpr int kNEq = 4;                     // HMOM
}

// C 链接入口, 给 Fluent 加载器
extern "C" {

DEFINE_INIT(mom_init, d) {
    // 从 case 读取物种、配置, 构建 thermo + model
    g_state = std::make_unique<ModelState>();
    // g_state->model = MOM::MakeAnyMomentMethod(thermo, "hmom");
    // 设置过程模型枚举...
}

DEFINE_ADJUST(mom_adjust, d) {
    if (!g_state) return;
    // 遍历所有单元, 每个算一次 ComputeCell, 源项存 C_UDMI
    thread_loop_c(t, d) {
        begin_c_loop(c, t) {
            real T = C_T(c,t), P = C_P(c,t);
            real mu = C_MU_L(c,t);
            // 取 Y_species, moments (从 UDS)
            // MOM::ComputeCell(model, T, P, Y, mu, moments);
            // auto src = MOM::GetSources(model);
            // for (i=0..NEq-1) C_UDMI(c,t,kUdmSootMomentSrcBase+i) = src[i];
        } end_c_loop(c,t)
    }
}

DEFINE_SOURCE(mom_source_m0, c, t, dS, eqn) {
    // 也可不缓存, 直接算; 这里读 ADJUST 预算好的值
    return C_UDMI(c,t,kUdmSootMomentSrcBase+0);
}
// ... mom_source_m1/m2/m3 同理

DEFINE_RW_FILE(mom_rw, fp) {
    // 读写模型状态用于 restart (可选)
}

} // extern "C"
```

### 步骤 4: 改写 `user_nt.udf`

```
# 选定 3ddp_node 文件夹下的 user_nt.udf
VERSION = 3ddp
PARALLEL_NODE = smpi
SOURCES = $(SRC)udf_bridge.cpp          # .cpp 扩展名 → cl 选 C++
USER_OBJECTS = $(LIBUDF)\MOM.lib          # 预编译静态库
# 需在 makefile 模板里加 /std:c++20 /EHsc /I<MOM include> /I<Eigen include>
```

注: 标准模板不直接暴露编译器选项行, 需深入 `makefile` 模板修改 `CSWITCH` / `CPPFLAGS`。实测需手动编辑 Fluent 安装目录下的 makefile 模板, 或用 `nmake` 命令行传 `CSWITCH`。

### 步骤 5: 处理并行

- `node_host_to_node_sync_*` 同步模型配置参数
- 每个计算节点独立 `g_state` (进程隔离, 天然安全)
- 共享内存多线程: 用 `thread_local` 或每线程实例 (Fluent 共享内存并行较少见, 通常 MPI 分区)

### 步骤 6: 注册 UDS 与源项

- Fluent 中定义 4 个 UDS (对应 HMOM 四矩)
- 每个 UDS 的源项指向对应 `DEFINE_SOURCE`
- 设定 UDS 对流/扩散项 (矩输运方程: 对流+扩散, 扩散系数可设 0 或小值)
- 初始化 UDS 值 (DEFINE_INIT 内)

### 步骤 7: 验证

- 用 `tests/test_verify_all_variants.cpp` 的已知算例做基准: 同一 T/P/Y 输入, Fluent 内调用 vs 离线调用, 源项数值应一致 (容差内)
- 先单机串行验证, 再上并行
- 检查质量/碳守恒 (源项积分 vs 碳烟总量变化)

---

## 6. 替代方案对比

| 方案 | 优点 | 缺点 |
|------|------|------|
| **A. 预编译 MOM.lib + C++ 桥接 (推荐)** | 性能最优, 全功能, 官方支持链接外部库 | 需改 makefile, ABI 对齐, 工程量大 |
| B. 把 MOM 源码直接加入 UDF SOURCES | 无链接问题 | 编译时间爆炸 (模板), 每次改 UDF 重编全部 |
| C. MOM 编译为 DLL, UDF 运行时 LoadLibrary | 解耦 | Fluent 加载 libudf.dll 时再 LoadLibrary 二级 DLL, 路径/依赖管理复杂 |
| D. 通过文件交换 (Fluent 写输入, 外部进程算, 读回) | 完全解耦 | 性能灾难 (I/O), 不实用 |
| E. 用 Fluent 的 Scheme/UDF 混合 + TCP socket | 跨语言灵活 | 延迟高, 复杂, 非标准 |

**推荐 A**: 平衡性能、功能与官方支持度。

---

## 7. 风险清单

| 风险 | 等级 | 缓解 |
|------|------|------|
| MSVC 版本/CRT 不匹配 → 链接失败或运行时崩溃 | 高 | 严格用 Fluent 捆绑 MSVC 重建 MOM |
| `udf.h` 在 C++20 下有兼容性警告/错误 | 中 | `-Wno-*` 或小改 udf.h 局部; 多数版本已兼容 |
| 并行下模型状态不同步 | 中 | DEFINE_ADJUST 内同步配置; 每节点独立实例 |
| UDS 数量上限 | 低 | Fluent 通常允许 ≥ 6 UDS, 够 HMOM(4) |
| 单位制不一致 (SI vs Fluent 默认) | 中 | 全程 SI, 源项转 kg/m³·s 等, 文档化换算 |
| Eigen 在 C++20 + MSVC 下的已知警告 | 低 | Eigen 3.4 兼容 C++20, 可能需抑制警告 |
| Fluent 版本升级改 makefile 模板 | 中 | 锁定版本, 升级时重验证 |

---

## 8. 何时值得做

**值得**:
- 长期项目, 需要在 Fluent 中反复跑碳烟/TiO₂ 模拟
- 已有 mom 验证算例, 想用 Fluent 的网格/湍流/辐射能力
- 团队有 C++ + CMake + Fluent UDF 经验

**不值得**:
- 一次性分析 — 直接用 OpenSMOKE++ 自家求解器或 mom 的测试驱动更省事
- 只需粗略源项 — Fluent 自带的碳烟模型 (Brookes-Moss 内置版) 即可
- 无 C++ 工具链经验 — 集成调试成本可能高于重写

---

## 9. 参考资料

- Ansys Fluent UDF Manual, DEFINE_SOURCE (node49.htm)
- Ansys Fluent UDF Manual, Build UDF Library (node133.htm)
- Ansys Fluent UDF Manual, Link Precompiled Object Files (node135.htm)
- mom 仓库: https://github.com/acuoci/mom.git
- mom 文档: 见 `DOCUMENTATION.md`
- 社区: cfd-online "Integration of a Custom C++ Model into FLUENT"; ResearchGate 静态库链接讨论

---

*评估基于官方文档与 mom 仓库架构分析; 实际可行性需在具体 Fluent 版本 + MSVC 工具链上原型验证。*
