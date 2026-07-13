# OpenSMOKE++ MOM — 架构与使用文档

> 仓库: https://github.com/acuoci/mom.git
> 作者: Alberto Cuoci, CRECK Modeling Lab, Politecnico di Milano
> 许可证: GPL-3.0-or-later
> 语言标准: C++20 (核心) / C++23 (字典解析可选)
> 构建系统: CMake (目标化, 4+ 兼容)

本文档是对 mom 仓库的系统性梳理，供后续二次开发、集成到 CFD 求解器、或扩展新矩方法变体时参考。

---

## 1. 项目定位

mom 是 OpenSMOKE++ 生态下的 **矩方法 (Method of Moments, MOM)** 库，用于在 CFD 计算中给出 **颗粒相 (soot / 金属氧化物) 的源项**。它不是一个独立求解器，而是一个被宿主 CFD 代码在每个网格单元调用、返回矩方程右端项 (源项 + 雅可比贡献) 的库。

典型调用链:
```
CFD 求解器 (每个单元/每步)
   ↓  传入 T, P, Y_species, mu, moments
MOM 库 (Compute / ComputeCell)
   ↓  返回 dM/dt 的源项 (std::span<const double>, NEq 维)
CFD 求解器组装到输运方程
```

支持四类颗粒系统:
| 变体 | NEq | 应用 | NDF 假设 | 关键过程 |
|------|-----|------|----------|----------|
| **HMOM** | 4 | 碳烟 | 双 δ 双峰 (核化小模态 + 大模态) | 成核/HACA生长/冷凝/氧化/凝并(9子机制) |
| **ThreeEquations** | 3 | 碳烟 | 对数正态 | 成核/HACA生长/氧化/凝并 |
| **BrookesMoss** | 2 | 碳烟 | (不重建NDF) | 成核/生长/氧化/凝并 |
| **MetalOxide** | 3 | TiO₂ | 对数正态 | 成核/凝并/**烧结**(Kruis) |

---

## 2. 目录结构

```
mom/
├── README.md, AI_CONTEXT.md, LICENSE, .clang-format
├── CMakeLists.txt
├── cmake/MOMConfig.cmake.in          # 安装配置模板
├── include/                          # 公共头 (对外可见)
│   ├── MOM/                          # 框架核心
│   │   ├── MOM.hpp                   # 总伞头
│   │   ├── MomentMethodConcept.hpp   # C++20 concept (权威契约)
│   │   ├── MomentMethodBase.hpp+.tpp # CRTP 基类
│   │   ├── AnyMomentMethod.hpp+.tpp  # std::variant 类型擦除
│   │   ├── MomVariantList.hpp        # 变体注册表 TypeList
│   │   ├── MomentMethodReporter.hpp  # 只观测者, 写输出列
│   │   ├── ProcessFlags.hpp          # 过程模型枚举
│   │   ├── ManifestReader.hpp        # 运行时清单 (manifest) 读取器 (仅声明, 见 §15)
│   │   └── ThermoProxy.hpp           # 热力学后端解耦 (BasicThermoData)
│   ├── HMOM/*.hpp+.tpp
│   ├── BrookesMoss/*.hpp+.tpp
│   ├── ThreeEquations/*.hpp+.tpp
│   └── MetalOxide/*.hpp+.tpp
├── src/                              # 私有实现 + 显式实例化
│   ├── {Variant}_BasicThermoData.cpp # extern template 实例化 (BasicThermoData)
│   ├── {Variant}_Grammar.cpp         # 仅 MOM_USE_DICTIONARY=ON 编译
│   ├── HMOM.cpp / BrookesMoss.cpp / ThreeEquations.cpp  # 非模板实现
│   ├── MOM/ThermoProxy_OpenSMOKE.cpp
│   └── Utilities/OutputFileColumns.cpp
└── tests/
    ├── test_verify_all_variants.cpp  # 物理正确性回归
    ├── test_loop_benchmark.cpp       # 性能基准
    ├── reporter_usage_example.cpp    # Reporter 用法示例
    ├── test_dictionary_setup.cpp     # 字典解析 (仅字典模式)
    ├── test_oh_fixture.cpp           # OH-only 氧化夹具 (G1 化学计量验证, 无 Fluent)
    └── g1_closure_analysis.py        # G1 闭合分析 (Python 3.13 stdlib, 只读)
```

**离线机制工具** (仓库根目录 `tools/`, 与 `mom/` 平级, 详见 §15):

```
tools/
├── hmom_case_schema.py               # HMOM 案例 TOML 配置的 pydantic v2 模型 (8 模型, extra='forbid')
├── prepare_hmom_case.py              # CLI: validate 子命令 — 校验案例配置
├── mechanism_qualifier.py            # CHEMKIN NASA-7 热力学解析器 + G0 资格判定 (qualify 子命令)
├── manifest_writer.py                # 清单写出器: 消费 config + G0 报告 → 确定性 JSON manifest (write 子命令)
├── examples/
│   └── hmom-case.toml                # 有效配置样例
└── tests/
    ├── test_config_schema.py         # 配置 schema 单元测试 (17)
    ├── test_mechanism_qualifier.py   # 机制资格判定单元测试 (22)
    ├── test_manifest_writer.py       # 清单写出器单元测试 (28)
    └── fixtures/sample_thermo.dat    # 9 物种 CHEMKIN 样例 (H,OH,O2,H2,H2O,CO,C2H2,A4,N2)
```

**头/源分离约定** (见 AI_CONTEXT.md):
- `include/` = 公共 API (`.hpp` 声明 + `.tpp` 模板实现)
- `src/` = 私有实现 + 显式模板实例化
- 通过 `extern template` 隐藏模板, 降低编译依赖与二进制膨胀
- 宏 `MOM_COMPILED_LIBRARY`: 库构建时定义, 头文件中 `.tpp` 由 `#if !defined(MOM_COMPILED_LIBRARY)` 守护, 避免重复实例化

---

## 3. 核心架构

### 3.1 CRTP 基类 — 零虚函数设计

```cpp
template <typename Derived, std::size_t NEq>
class MomentMethodBase {
protected:
    using MomentVector = Eigen::Matrix<double, NEq, 1>;  // 定长, 栈上, SIMD
    MomentVector moments_;                                // 当前单元矩
    Eigen::VectorXd omega_gas_;                           // 气相源项 (setup 时定容)
public:
    // 每个过程源项: if constexpr 检测派生类是否实现 _impl()
    auto sources_nucleation() const {
        if constexpr (requires(const Derived& d){ d.sources_nucleation_impl(); })
            return static_cast<const Derived*>(this)->sources_nucleation_impl();
        else
            return std::span<const double>{kZeroData, NEq};  // 静态 .rodata 零数组
    }
    // ... coagulation / condensation / growth / oxidation / sintering 同构
};
```

要点:
- **零虚函数** → 编译期多态, 派生类方法 `[[gnu::always_inline]]` 即使 `-O0` 也内联
- **`if constexpr (requires(...))` 检测** → 派生类不实现某过程时, 基类返回指向静态零数组的 `std::span` (零拷贝, 零分配)
- **定长 Eigen 向量** → 栈上, 无堆分配, 友好 SIMD 向量化
- **`std::span<const double>` 返回** → 调用方零拷贝读取源项

### 3.2 C++20 Concept — 权威公共契约

`MomentMethodConcept.hpp` 定义 `MomentMethod<M>` concept, 是判断一个类型是否为合法矩方法的**唯一权威**。另有 `HasReconstructedNDF<M>` concept (要求 `ReconstructedNDF` + `ReconstructedNormalizedNDF` 成员)。

自由函数 `ComputeCell(model, T, P, Y, mu, moments) noexcept` 是单单元调用的便捷入口。

### 3.3 热力学后端解耦 — ThermoMap concept

`ThermoProxy.hpp` 定义 `ThermoMap<T>` concept, 使矩方法与具体热力学库解耦。库自带满足该 concept 的 `BasicThermoData` 结构体:
- 物种名列表、分子量 (kg/kmol)
- 原子计数 `nc/nh/no/nn/nti` (碳/氢/氧/氮/钛)
- 索引 0-based, `IndexOfSpecies` 缺失返回 -1

`ThermoProxy_OpenSMOKE.cpp` 提供 OpenSMOKE++ 后端适配 (仅集成 OpenSMOKE++ 时使用)。**新后端只需满足 ThermoMap concept 即可接入, 无需改动矩方法代码。**

---

## 4. 四个变体详解

所有变体在 `MomVariantList.hpp` 中通过 `detail::TypeList<HMOM,BrookesMoss,ThreeEquations,MetalOxide>` 注册, 供 `AnyMomentMethod` 的 `std::variant` 枚举。

### 4.1 HMOM (Hybrid Method of Moments) — NEq=4

- **矩向量**: `[M00, M10, M01, N0]` (归一化, mol/m³)
- **NDF**: 双 δ 双峰 — 小模态 (核化节点 `N0@V0`), 大模态 (`NL@VL`)
- **建模过程**: 成核 + 表面生长 (HACA 5步) + 冷凝 + 氧化 + 凝并
- **凝并 9 子机制**: 离散 SS/SL/LL + 连续 SS/SL/LL + 总计 (用于诊断分解)
- **不建模**: 烧结
- **标签**: `{"HMOM","hmom"}`
- **配置枚举**: `StickingModel`, `FractalDiameterModel`, `CollisionDiameterModel`
- **NDF 重建**: `NDFReconstructionData`, `kNDFSmearSigmaLnV=0.5` (对数正态涂抹, 仅可视化)
- **HACA Arrhenius 参数**: `A1f..A5, eff6`
- **氧化化学计量 (per-C2)**: 氧化速率 `R_oxid_C2` 为 **每 C2 对** 速率 (每单位移除 2 个碳原子, 由 `S_M10` 碳损失独立验证)。两条氧化通道均按 per-C2 配平:
  - O2 通道: `C2 + O2 → 2CO` (1 O2, 2 CO)
  - OH 通道: `C2 + 2OH → 2CO + H2` (2 OH, 2 CO, 1 H2)
  - **G1 修复 (HMOM.tpp:1103-1124)**: OH 通道此前误用 per-C 系数 (1 OH, 1 CO, 0.5 H2) 施加于 per-C2 速率, 导致 OH/CO/H2 各欠计 2×、质量与碳约 50% 不闭合、放热低估 50%。修复后气相项乘 2 (OH/CO) 与 0.5→1.0 (H2), 与 O2 通道一致; 质量/C/H/O/焓闭合均 < 1e-15。仅改化学计量, 不动 `kox_OH_` 速率常数前因子。验证见 `.swarm/evidence/1/`。
- **关键私有方法**: `DimerConcentration`, `SootKineticConstants`, `SootNucleationM4`, `SootSurfaceGrowthM4`, `SootOxidationM4`, `SootCondensationM4`, `SootCoagulationM4` (+ SmallSmall/SmallLarge/LargeLarge/Continuous 变体), `GetMissingMoment`, `GetBetaC`, `PAHDimerizationRate`
- **Reporter 扩展钩子**:
  - `variant_prefix_output`: np, ss, vs, N0, NL, alphaL, dpL, npL, d63, sigma_dp/np, sigma_dp/np_L, gsd_*, omega*
  - `variant_suffix_output`: 9×NEq 凝并子向量 (ScoaTot/Dis/DisSS/DisSL/DisLL/Con/ConSS/ConSL/ConLL)
  - `ndf_extra_output`: N0, V0, dp0, NL, VL, dpL_mean, sigma_ndf
- **文献**: Mueller 2009; Attili & Bisetti 2013/2014

### 4.2 ThreeEquations — NEq=3

- **矩向量**: `[Ys, NsNorm, Ss]`
- **NDF**: 对数正态, 实现 `HasReconstructedNDF`
- **建模**: HACA 成核/生长/氧化/凝并
- **文献**: Eberle 2015

### 4.3 BrookesMoss — NEq=2

- **矩向量**: `[Ys, bs]`
- **最简变体**, 无冷凝、无烧结 (零跨度返回)
- **不重建 NDF** (无 `HasReconstructedNDF`)
- **文献**: Brookes & Moss 1999

### 4.4 MetalOxide (TiO₂) — NEq=3

- **矩向量**: `[YTiO2, NTiO2N, STiO2]`
- **建模**: 成核 + 凝并 + **烧结 (Kruis 粘性流模型)**
- **不建模**: 表面生长、氧化
- `PlanckCoeffModel::None` (介电颗粒, 无辐射吸收)
- `rho_particle` 默认 3900 kg/m³
- 实现 `HasReconstructedNDF`
- **文献**: Pratsinis 1988; Jossen 2005

---

## 5. 类型擦除与运行时选择

### 5.1 AnyMomentMethod

```cpp
template <typename Thermo>
using AnyMomentMethod = AllVariants::AsVariant<Thermo>;
// = std::variant<HMOM<Thermo>, BrookesMoss<Thermo>, ThreeEquations<Thermo>, MetalOxide<Thermo>>
```

工厂函数:
```cpp
auto model = MakeAnyMomentMethod(thermo, "hmom");  // 标签匹配, 递归 FactoryHelper
```
匹配每个变体的静态 `variant_labels` (大小写不敏感, 见 HMOM 的 `{"HMOM","hmom"}`)。

### 5.2 类型擦除辅助函数

| 函数 | 说明 |
|------|------|
| `GetNEquations(model)` | 返回当前变体 NEq (一次 `std::visit`) |
| `SetState(model, T, P, Y)` | 设置单元热力学状态 |
| `SetViscosity(model, mu)` | 设置气相粘度 |
| `SetMoments(model, moments)` | 设置当前矩 |
| `Compute(model)` | 触发源项计算 |
| `ComputeCell(model, T, P, Y, mu, moments)` | 单单元一站式 (`std::visit` 一次) |
| `ForEachCell(model, callback)` | **循环外做 `std::visit`**, 回调内复用活动变体指针 — 关键性能优化 |
| `GetSources(model)` / `GetOmegaGas(model)` / `GetVolumeFraction(model)` | 结果取值 |

**性能要点**: 逐单元调用时优先用 `ForEachCell`, 把 `std::visit` 的分支预测开销移到循环外; 避免在 cell 循环内反复 `std::visit`。

### 5.3 字典解析 (可选, C++23)

仅 `MOM_USE_DICTIONARY=ON` 时编译 `*_Grammar.cpp`, 链接 `libdictionary_parser.a`。使用 C++23 `std::expected` 做错误处理。`SetupFromDictionary` 从 OpenSMOKE++ 字典输入配置变体。

---

## 6. Reporter — 只读输出列观测者

`MomentMethodReporter.hpp` 是与求解器解耦的输出层, 负责按列写出诊断量。

### 6.1 构造与使用

```cpp
MomentMethodReporter reporter{output_columns, species_names};
reporter.WriteHeader(model);   // 模板上, 一次表头
// ... 每步:
reporter.WriteRow(model);      // 写一行数据
reporter.WriteReconstructedNDF(model);  // 仅 HasReconstructedNDF, 否则静默 no-op
```

对 `AnyMomentMethod` 有重载 (每次调用一次 `std::visit`, 非 per-cell)。

### 6.2 列布局 (顺序)

```
Block1 核心: Ys, Ns, Ss, fv, dp, dc
   → variant_prefix_output 钩子 (变体自定义列)
Block3 输运: D
Block4: Sall(j)  (总源项各分量)
Block5 分过程: Snuc/Scoa/Scon/Sgro/Soxi/Ssin
   (未拥有的过程标 [ZF] 标签, 通过 requires _impl 检测)
   → variant_suffix_output 钩子
```

### 6.3 变体扩展协议

变体实现以下模板方法即可注入自定义列, Reporter 通过 `if constexpr requires` 自动检测:
```cpp
template <typename CB>
void variant_prefix_output(CB&& cb) const;  // cb(string_view label, double value)
template <typename CB>
void variant_suffix_output(CB&& cb) const;
template <typename CB>
void ndf_extra_output(CB&& cb) const;
```
**同一实现同时服务于 header(label 模式) 和 row(value 模式)**, 避免重复代码。

### 6.4 重建 NDF 输出

`WriteReconstructedNDF` 要求 `HasReconstructedNDF` (BrookesMoss 静默跳过):
- 6 核心列: `nu, n, nbar, dsph, nd, ndbar` (单位 nm)
- + `ndf_extra_output` 钩子
- 雅可比 `|dv/dd_sph| = (π/2)·d²`

---

## 7. 过程模型枚举 (ProcessFlags.hpp)

全部 `enum class : int`, 通用取值 `Off=0 / Standard=1 / Extended=2` (适用时):

| 枚举 | 取值 |
|------|------|
| `NucleationModel` | Off / Standard / Extended |
| `CoagulationModel` | Off / Standard / Extended |
| `SurfaceGrowthModel` | Off / Standard / Extended |
| `OxidationModel` | Off / Standard / Extended |
| `CondensationModel` | Off / Standard / Extended |
| `SinteringModel` | Off / Standard / Extended |
| `ThermophoreticModel` | Off / Standard / Extended |
| `PlanckCoeffModel` | **None / Smooke / Kent / Sazhin** |

字符串解析器: `PlanckCoeffModelFromString`, `NucleationModelFromString`, `OxidationModelFromString`。

---

## 8. 物理常数与辐射

### 8.1 物理常数 (MomentMethodBase, CODATA 2018 / IUPAC 2021)

| 符号 | 值 | 单位 |
|------|----|------|
| `kB_` | 1.380649e-23 | J/K |
| `Nav_mol_` | 6.02214076e23 | 1/mol |
| `Nav_kmol_` | 6.02214076e26 | 1/kmol |
| `Rgas_` | 8314.46261815324 | J/(kmol·K) |
| `WC_` | 12.011 | kg/kmol (碳) |
| `pi_`, `sqrt2_` | — | — |

### 8.2 Planck 辐射吸收系数 (`planck_coefficient()`, switch 选择)

| 模型 | 公式 | 文献 |
|------|------|------|
| Smooke | `1232.4 · fv · T` | Smooke 2005 |
| Kent | `1.3e5 · fv · T` | Kent & Honnery 1990 |
| Sazhin | `fv · (630 + 0.63·T − 1e-4·T²)` | Sazhin 1994 |
| None | 0 | (介电颗粒, 如 TiO₂) |

`fv` = 体积分数, `T` = 气相温度 [K]。

---

## 9. CMake 构建系统

### 9.1 选项

| 选项 | 默认 | 说明 |
|------|------|------|
| `MOM_BUILD_SHARED` | OFF | 构建共享库 (否则静态) |
| `MOM_BUILD_TESTS` | ON | 构建 tests/ |
| `MOM_INSTALL` | ON | 生成 install 目标 |
| `MOM_USE_DICTIONARY` | OFF | 启用字典解析 (需 `DICTIONARY_DIR`, 升级到 C++23) |

### 9.2 目标

- **`MOM_headers`** (INTERFACE): 纯头文件使用方式, 仅 `cxx_std_20` + include 目录 + Eigen 链接
- **`MOM`** (STATIC/SHARED): 编译库, 定义 `MOM_COMPILED_LIBRARY` → 头文件启用 `extern template` 抑制重复实例化

`mom_configure_target()` 统一设置: include 目录 (`BUILD_INTERFACE`/`INSTALL_INTERFACE`) + 链接 `Eigen3::Eigen` + `cxx_std_20`。

### 9.3 Eigen 依赖解析顺序

1. `EIGEN3_DIR` 显式覆盖 → `find_package`
2. 否则 `find_package(Eigen3)`
3. 否则 `FetchContent` 从 gitlab 拉取 3.4.0

### 9.4 字典模式特殊处理

`MOM_USE_DICTIONARY=ON` + `DICTIONARY_DIR` 设置时:
- 编译 `*_Grammar.cpp` (拉 OpenSMOKE++ 头)
- 链接 `${DICTIONARY_DIR}/lib/libdictionary_parser.a`
- include `${DICTIONARY_DIR}/include`
- 定义 `MOM_USE_DICTIONARY`, **升级到 `cxx_std_23`** (用 `std::expected`)

### 9.5 源文件清单 (始终编译)

- `{ThreeEquations,HMOM,BrookesMoss,MetalOxide}_BasicThermoData.cpp` (显式模板实例化)
- `Utilities/OutputFileColumns.cpp`
- `HMOM.cpp / BrookesMoss.cpp / ThreeEquations.cpp` (非模板实现)

### 9.6 安装

GNUInstallDirs + `MOMTargets` 导出 (命名空间 `MOM::`) + `cmake/MOMConfig.cmake.in` 模板, 支持下游 `find_package(MOM)`。

### 9.7 编译警告

GCC/Clang: `-Wall -Wextra -Wpedantic`。`CMAKE_CXX_EXTENSIONS OFF` (禁用编译器扩展)。

### 9.8 源/执行字符集 (MSVC `/utf-8`)

MOM 源文件为 UTF-8 编码 (注释/字符串含 `C₂H₂`、`→`、`✓` 等非 ASCII 字符)。MSVC 默认按宿主代码页 (如 936/GBK) 解码源文件, 会破坏这些多字节序列并触发伪 C2065/C3688 错误。`CMakeLists.txt` 通过 `add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/utf-8>)` 强制源与执行字符集均为 UTF-8 (仅对 MSVC 生效)。

---

## 10. 代码规范 (AI_CONTEXT.md)

| 规则 | 说明 |
|------|------|
| 标准 | C++20 核心, C++23 仅字典模式 |
| 复杂模板函数 | 尾置返回类型 (`auto f() -> T`) |
| `[[nodiscard]]` | 非 void 无副作用函数 |
| `noexcept` | 强异常安全 (move/getter/数学计算) |
| `constexpr` | 尽可能激进使用 |
| 头/源分离 | `include/` 公共, `src/` 私有 + 显式实例化 |
| `extern template` | 隐藏模板, 降低编译耦合 |
| 依赖最小化 | 核心仅 stdlib + Eigen; 重生态隔离 |
| 缓存友好 | `std::vector`, 定长 Eigen 栈向量 |
| 只读字符串 | `std::string_view` |
| 移动语义 | 优先 |
| CMake | 目标化 (`target_*`), 禁遗留全局变量 |

### .clang-format

LLVM 基础, Allman 大括号, 4 空格缩进, 100 列, 指针左对齐, `SortImports: Never`, 最低 clang-format 14。
Include 顺序: STL → Eigen → 项目内。

---

## 11. 典型集成示例

### 11.1 静态多态 (编译期已知变体, 最高性能)

```cpp
#include "MOM/MOM.hpp"

MOM::HMOM<MOM::BasicThermoData> model(thermo_data);
model.set_nucleation_model(MOM::NucleationModel::Standard);
model.set_planck_coefficient_model(MOM::PlanckCoeffModel::Smooke);

for (each cell) {
    MOM::ComputeCell(model, T, P, Y_gas, mu, moments_span);
    auto src = model.sources_total();   // std::span<const double, 4>
    // 组装到输运方程
}
```

### 11.2 运行时选择 (配置文件决定变体)

```cpp
auto model = MOM::MakeAnyMomentMethod(thermo, "hmom");  // std::variant

MOM::ForEachCell(model, [&](auto& concrete) {
    // 此处 concrete 已是具体变体类型, 无额外 std::visit 开销
    MOM::SetState(concrete, T, P, Y);
    MOM::SetViscosity(concrete, mu);
    MOM::SetMoments(concrete, moments);
    MOM::Compute(concrete);
    auto src = MOM::GetSources(concrete);
    // ...
});
```

### 11.3 输出诊断

```cpp
MOM::OutputFileColumns cols;
MOM::MomentMethodReporter reporter{cols, species_names};
reporter.WriteHeader(model);
// 每步:
reporter.WriteRow(model);
reporter.WriteReconstructedNDF(model);
```

---

## 12. 扩展指南: 添加新变体

1. 在 `include/NewVariant/` 新建 `.hpp` + `.tpp`, 定义:
   ```cpp
   template <typename Thermo>
   class NewVariant : public MOM::MomentMethodBase<NewVariant<Thermo>, NEq> {
       static constexpr std::array variant_labels{"NewVariant","newvariant"};
       // 实现 sources_X_impl() for owned processes
       // 可选: variant_prefix_output / variant_suffix_output / ndf_extra_output
       // 可选: ReconstructedNDF / ReconstructedNormalizedNDF (满足 HasReconstructedNDF)
   };
   ```
2. 在 `MomVariantList.hpp` 的 `TypeList` 末尾追加 `NewVariant`
3. 在 `src/` 添加 `NewVariant_BasicThermoData.cpp` (显式实例化 `NewVariant<BasicThermoData>`)
4. 在 `CMakeLists.txt` 源文件清单加入新 `.cpp`
5. (可选) 添加 `NewVariant_Grammar.cpp` (字典模式)
6. 在 `tests/test_verify_all_variants.cpp` 加入回归用例

**无需改动 Reporter / AnyMomentMethod / Concept** — `if constexpr requires` 与 `TypeList` 自动接入。

---

## 13. 参考文献

| 主题 | 文献 |
|------|------|
| HMOM | Mueller, 2009; Attili & Bisetti, 2013/2014 |
| ThreeEquations | Eberle, 2015 |
| BrookesMoss | Brookes & Moss, 1999 |
| TiO₂ | Pratsinis, 1988; Jossen, 2005 |
| Planck 辐射 | Smooke 2005; Kent & Honnery 1990; Sazhin 1994 |

---

## 14. 未深入阅读的部分 (后续按需查阅)

- `.tpp` 实现体 (MomentMethodBase.tpp, HMOM.tpp, AnyMomentMethod.tpp, 各变体 .tpp) — 具体物理公式
- `src/*.cpp` 物理实现 (HACA 动力学、凝并核、烧结模型等)
- `tests/*.cpp` — 回归基准的具体测试用例与期望值
- `OutputFileColumns.h` — 列管理的完整接口
- BrookesMoss / ThreeEquations / MetalOxide 的 `.hpp` 完整内容 (已确认 NEq 与 concept 契合, 细节未逐行读)
- `*_Grammar.cpp` — OpenSMOKE++ 字典语法

本文档基于对仓库结构、公共头、CMake、README、AI_CONTEXT 的通读整理, 足以支撑集成与二次开发决策; 物理公式的精确形式需查阅对应 `.tpp`/`.cpp` 与原始文献。

---

## 15. 离线机制工具与运行时清单 (Offline Tools & Manifest)

`tools/` (仓库根目录, 与 `mom/` 平级) 存放 **每个案例运行前执行一次** 的 Python 离线工具链。它把人工编写的案例配置和 CHEMKIN 机制转换为一份 **确定性、带完整性校验的运行时清单 (manifest)** — 一个 JSON 文件, 由 C++ UDF 桥接层 (`MOM::ManifestReader`) 在求解器启动时读取。

工具链为 **Python 3.11+**, 除配置校验器与清单写出器用到 **pydantic v2** 外, 全部仅依赖标准库。每个工具既可作脚本 (`python tools/<tool>.py …`) 也可作模块 (`python -m tools.<tool> …`) 运行。

### 15.1 三阶段流水线

```
hmom-case.toml ──▶ prepare_hmom_case.py validate   (schema / 语义校验)
                          │
thermo.dat ──────▶ mechanism_qualifier.py qualify  (CHEMKIN NASA-7 解析 + G0 闸门)
                          │                              └──▶ g0_report.json
                          ▼
        manifest_writer.py write config.toml g0_report.json  (交叉校验 + 写出)
                          │
                          ▼
                  manifest.json  ──▶  MOM::ManifestReader::load()  (C++ 求解器启动)
```

### 15.2 案例配置校验器 — `prepare_hmom_case.py`

校验 UTF-8 TOML 案例配置, 使用严格 pydantic v2 schema (`hmom_case_schema.py`, 8 个模型, 全部 `extra="forbid"` — 拒绝未知键、类型不符、缺失字段)。schema 强制 v1 不变量: `schema_version = 1`、`bulk_species = "N2"`、`thermophoretic_model = 0`、`uds_count = 4`、64 位十六进制 `expected_mechanism_sha256`、`pah_species` 不得与保留气相物种名冲突。

配置分 8 个必需节: `[case] [mechanism] [hmom] [fluent] [radiation] [energy] [validation]` + 顶层 `schema_version`。`[fluent]` 节由 `udm_mode` 派生 `udm_count` — release 模式 **18** 槽, debug 模式 **38** 槽。

```bash
python tools/prepare_hmom_case.py validate tools/examples/hmom-case.toml
```

`validate` 子命令: 成功打印摘要并退出 0; 失败按字段将错误打印到 stderr 并退出 1。样例见 `tools/examples/hmom-case.toml`。

### 15.3 机制资格判定 (G0) — `mechanism_qualifier.py`

解析 CHEMKIN `THERMO` 段 (4 行 × 80 列 NASA-7 条目), 产出 **G0 资格报告**。G0 检查:

- 固定必需物种齐备 — `H, OH, O2, H2, H2O, CO, C2H2, N2` — 外加可配置 PAH 物种;
- `N2` 是机制中 **最后** 一个物种;
- PAH 物种索引独立, 且不与保留名冲突;
- 全部 NASA-7 系数有限、分子量为正、生成焓有限;
- 无重复物种名。

对每个物种由 NASA-7 多项式计算分子量、标准生成焓 `Hf0(298.15 K)`、比热 `Cp(298.15 K)`, 并记录温度区间与原始系数。

```bash
python tools/mechanism_qualifier.py qualify tools/tests/fixtures/sample_thermo.dat \
    --pah-species A4 --output g0_report.json
```

选项: `--pah-species` (默认 `A4`)、`--output <file>` (写 JSON 报告, 否则打印到 stdout)、`--quiet` (抑制 JSON 输出)。合格退出 0, 否则退出 1。报告含 `thermo_sha256`、`species_order`、逐物种记录, 以及解释失败原因的 `gaps` 列表。

### 15.4 清单写出器 — `manifest_writer.py`

同时消费已校验的 config 与 G0 报告, **交叉校验** 二者 (匹配 `pah_species`、`thermo_file`、`N2` 末位、必需物种齐备, 并要求 `qualified = true`), 写出 **确定性** JSON 清单。确定性由 `sort_keys=True` 序列化 + 内嵌自 SHA-256 (在 `manifest_sha256` 字段清空后对清单计算) 保证。

```bash
python tools/manifest_writer.py write tools/examples/hmom-case.toml g0_report.json \
    --output manifest.json
```

`write` 子命令: 成功退出 0。带 `--output` 时打印物种数与 `manifest_sha256`; 否则写到 stdout。

**样例夹具**: `tools/tests/fixtures/sample_thermo.dat` 是 9 物种 CHEMKIN 样例 (`H, OH, O2, H2, H2O, CO, C2H2, A4, N2`), 以 `--pah-species A4` 通过 G0, 便于无生产机制时演练完整流水线。单元测试位于 `tools/tests/` (`test_config_schema.py` 17, `test_mechanism_qualifier.py` 22, `test_manifest_writer.py` 28)。

### 15.5 清单 JSON 格式

`manifest_writer.py` 产出的清单为顶层对象, 关键字段:

| 字段 | 类型 | 说明 |
|------|------|------|
| `magic` | string | 固定 `"HMOM-MANIFEST"` |
| `schema_version` | int | 固定 `1` |
| `generator_version` | string | 写出器版本 (如 `"1.0.0"`) |
| `generated_at` | string | ISO 8601 UTC 时间戳 (仅溯源, **不参与** 确定性契约) |
| `sources` | object | `config_file / config_sha256 / thermo_file / thermo_sha256` |
| `species_count` | int | 物种条目数 |
| `species` | array | 物种块 (机制顺序), 每项见下 |
| `indices` | object | `pah_index`、`bulk_index`、`source_species` (映射 C2H2/H2/O2/OH/CO/PAH → 机制索引) |
| `hmom_options` | object | 解析后的 HMOM 过程/参数 (16 项) |
| `fluent_options` | object | `udm_mode / udm_count / uds_count` |
| `udm_layout` | object | 固定 UDM 槽布局 (见 §15.6) |
| `operating_pressure_pa` | float | 固定操作压力 [Pa] |
| `radiation` | object | `enabled / model / zones / require_gray_wsggm_do / require_native_soot_off` |
| `energy_options` | object | `reaction_heat_enabled / soot_formation_enthalpy_j_kg / require_sensible_enthalpy_equation` |
| `manifest_sha256` | string | 对清单 (该字段清空后) 的 SHA-256, 完整性/确定性保证 |

**物种块 (`species[]`)** 每项字段: `name`、`mw_kg_kmol`、`mw_kg_mol`、原子计数 `nc/nh/no/nn`、生成焓 `h_f0_j_kg` 与 `h_f0_j_mol`、`nasa7_low[7]`、`nasa7_high[7]`、温度界 `temp_low/temp_common/temp_high`。

**自 SHA-256 契约**: 将 `manifest_sha256` 置空 → 以 `sort_keys=True, indent=2` 序列化 → 对 UTF-8 字节取 SHA-256。`verify_self_sha256()` 据此校验一致性。

### 15.6 UDM 槽布局 (release 模式)

清单 `udm_layout` 记录 Fluent 用户自定义存储 (UDM) 的固定槽位 (release 模式 18 槽; debug 模式 38 槽):

| 槽位 | 内容 |
|------|------|
| 0–3 | 矩源项 `moment_sources` (对应 HMOM NEq=4) |
| 4–9 | 物种源项 `species_sources`: C2H2=4, H2=5, O2=6, OH=7, CO=8, PAH=9 |
| 10 | `bulk_dummy` (N2 占位) |
| 11 | `fv` — 体积分数 |
| 12 | `dp` — 一次颗粒直径 |
| 13 | `dc` — 碰撞直径 |
| 14 | `gamma_mom` |
| 15 | `q_hmom` — HMOM 放热 |
| 16 | `cache_epoch` |
| 17 | `status_bits` |

### 15.7 ManifestReader.hpp — C++ 清单读取器 (仅声明)

`include/MOM/ManifestReader.hpp` 是 C++20 纯声明头 (实现在 task 3.1 提供), 供 Fluent UDF 桥接层在启动时加载并校验清单。三个数据结构 + 一个读取器:

- **`struct ManifestSpecies`**: 单物种记录, 镜像清单 `species` 数组一项。字段: `name`、`mw_kg_kmol`、`mw_kg_mol`、`nc/nh/no/nn`、`h_f0_j_kg`、`h_f0_j_mol`、`std::array<double,7> nasa7_low`、`nasa7_high`、`temp_low/temp_common/temp_high`。
- **`struct Manifest`**: 完整解析后的清单。含头/溯源 (`magic`、`schema_version`、`generator_version`、`generated_at`)、来源 (`config_file/sha256`、`thermo_file/sha256`)、物种表 (`species_count` + `std::vector<ManifestSpecies> species`)、索引映射 (`pah_index`、`bulk_index`、`c2h2_index`、`h2_index`、`o2_index`、`oh_index`、`co_index`)、16 项 HMOM 选项、Fluent 选项 (`udm_mode`、`udm_count`、`uds_count`)、`operating_pressure_pa`、辐射 (`radiation_zones`、`radiation_enabled`、`radiation_model`)、能量选项, 以及 `manifest_sha256`。
- **`class ManifestParseError : public std::runtime_error`**: 解析/校验失败时抛出。
- **`class ManifestReader`**: 静态接口 —
  - `static Manifest load(const std::string& path)` — 从文件加载并校验;
  - `static Manifest parse(const std::string& json_text)` — 从内存 JSON 串解析并校验;
  - `static void validate(const Manifest& m)` — 校验结构不变量 (magic、schema 版本、自 SHA-256、物种数、N2 末位、源物种索引互异、`udm_count` 与 `udm_mode` 一致);
  - 常量: `kMagic = "HMOM-MANIFEST"`、`kSchemaVersion = 1`、`kUdsCount = 4`、`kUdmReleaseCount = 18`、`kUdmDebugCount = 38`。

---

## 16. C++ 上下文与适配层 (Context & Adapter Layer)

上一节的离线工具产出确定性 `manifest.json`。C++ **上下文与适配层** 在求解器启动时消费该清单, 将其转化为一个即用的 MOM 模型。它是 Python 产出的清单与编译期 MOM 内核之间的桥梁, 也是 Fluent UDF (或任意宿主 CFD 求解器) 在初始化阶段按 zone 调用一次的代码。

该层由三个可组合部件构成, 每个都配有对应的单元测试程序:

```
manifest.json
     │  MOM::ManifestReader::load()          (严格 JSON 解析 + fail-fast 校验)
     ▼
MOM::Manifest  (已解析、已校验的结构体)
     │  MOM::MechanismAdapter                 (Manifest → thermo + config + descriptor)
     ▼
BasicThermoData + HMOM::Config + MechanismDescriptor
     │  MOM::MomContext(manifest_path)        (拥有 thermo + HMOM 模型 + scratch 缓冲)
     ▼
MomContext::compute(T, P, Y, mu, moments)     (noexcept 逐单元入口)
```

### 16.1 ManifestReader — 严格 JSON 解析器

`include/MOM/ManifestReader.hpp` / `src/MOM/ManifestReader.cpp` 使用内置的 **nlohmann/json** 库解析由 `tools/manifest_writer.py` 写出的清单 (task 2.3 声明, 现已实现)。解析 **严格且 fail-fast**: 任何结构或语义违规均抛出 `ManifestParseError` (继承自 `std::runtime_error`) 并附带描述性消息, 而非返回部分填充的结构体。

```cpp
#include "MOM/ManifestReader.hpp"

MOM::Manifest m = MOM::ManifestReader::load("manifest.json");   // 任何违规即抛出
// 或从内存字符串:
MOM::Manifest m = MOM::ManifestReader::parse(json_text);
```

读取器在返回前强制 **11 项校验不变量**:

| # | 不变量 | 拒绝原因 |
|---|--------|----------|
| 1 | `magic == "HMOM-MANIFEST"` | 非 HMOM 清单 |
| 2 | `schema_version == 1` | 不支持的 schema |
| 3 | `species_count == species.size()` | 计数/数组不符 |
| 4 | `species_count > 0` | 空机制 |
| 5 | 末位物种为 `N2` | bulk 物种非末位 |
| 6 | 每个 `mw_kg_kmol > 0` | 非物理分子量 |
| 7 | 每个 `h_f0_j_kg` 有限 | 生成焓 NaN/Inf |
| 8 | 所有 `nasa7_low[7]` / `nasa7_high[7]` 有限 | 多项式系数 NaN/Inf |
| 9 | `uds_count == 4` | 输运标量数错误 |
| 10 | `udm_count` 与 `udm_mode` 一致 (release 18 / debug 38) | UDM 布局不符 |
| 11 | `operating_pressure_pa > 0` | 非物理操作压力 |

对宿主代码暴露静态契约常量: `kMagic`、`kSchemaVersion = 1`、`kUdsCount = 4`、`kUdmReleaseCount = 18`、`kUdmDebugCount = 38`。

### 16.2 MechanismAdapter — Manifest → 模型输入

`include/MOM/MechanismAdapter.hpp` / `src/MOM/MechanismAdapter.cpp` 将已校验的 `Manifest` 翻译为 HMOM 内核所需的三个对象:

- **`BasicThermoData`** — 物种名、**换算为 kg/kmol** 的分子量、每物种原子计数 (`nc/nh/no/nn`), 满足 `ThermoMap` concept。
- **`HMOM::Config`** — 由 `hmom_options` 映射的全部 **16 项** HMOM 过程/参数选项, 其中 `pah_species` 由清单 `pah_index` 解析得到。
- **`MechanismDescriptor`** — 上述两者未持有的其余物理数据: 每物种生成焓 `h_f^0`、NASA-7 低/高温多项式系数、源物种索引, 以及操作压力。

```cpp
MOM::MechanismAdapter adapter(manifest);
MOM::BasicThermoData     thermo     = adapter.thermo();
MOM::HMOM<...>::Config   config     = adapter.config();
MOM::MechanismDescriptor descriptor = adapter.descriptor();
```

**护栏 (Guard rails)**。适配器对 C++ 内核不支持的配置拒绝构建模型 (抛出而非静默降级):

- 拒绝 `sticking_model == "golaut"` (不支持的粘附模型);
- 拒绝 `thermophoretic_model != 0` (本桥接层不建模热泳)。

### 16.3 MomContext — 拥有型运行时上下文

`include/MOM/MomContext.hpp` / `src/MOM/MomContext.cpp` 提供宿主求解器在案例生命周期内持有的单一对象。`MomContext` **拥有** 热力学数据与 HMOM 模型, 外加可复用的 scratch 缓冲:

```cpp
struct MomContext {
    BasicThermoData        thermo_;   // 首先声明 —— 生命周期长于 model_
    std::unique_ptr<HMOM<BasicThermoData>> model_;
    // + MechanismDescriptor、scratch 缓冲、manifest hash、状态
};
```

**声明顺序是刻意的**: `thermo_` 在 `model_` 之前声明, 由于 HMOM 模型持有对热力学后端的引用, 因此析构时热力学数据保证晚于模型销毁。`MomContext` **不可拷贝、不可移动** (拷贝与移动均 `= delete`) —— 它是一个稳定、钉住的所有者。

构造函数接收清单路径并执行完整流水线 (加载 → 校验 → 适配 → 构造 HMOM)。逐单元入口为 `noexcept`:

```cpp
MOM::MomContext ctx("manifest.json");        // 加载 + 构建; 设置 is_valid()/error_message()
if (!ctx.is_valid()) { /* 处理 ctx.error_message() */ }

ctx.compute(T, P, Y.data(), mu, moments_span);   // noexcept —— 包装 MOM::ComputeCell
auto src = ctx.sources();                         // 零拷贝 span, 指向模型内部
```

`compute()` 按正确顺序设置状态/矩/粘度并调用 `MOM::ComputeCell`, 绝不跨越 C 边界抛出异常。上下文暴露 **12 个访问器**:

| 访问器 | 返回 |
|--------|------|
| `sources()` | 矩源项总向量 (span) |
| `omega_gas()` | 气相物种源项 |
| `V0()` | 核化节点颗粒体积 |
| `particle_density()` | 碳烟材料密度 [kg/m³] |
| `diffusion_coefficient()` | 有效扩散系数 (> 0) |
| `descriptor()` | `MechanismDescriptor` |
| `zone_ids()` | 清单中的 Fluent zone ID |
| `manifest_hash()` | 已加载清单的 `manifest_sha256` |
| `is_valid()` | 构造是否成功 |
| `error_message()` | 失败描述 (有效时为空) |
| `n_species()` | 物种数 |
| `species_names()` | 物种名列表 |

### 16.4 单元测试

三个测试程序 (**共 28 个测试用例**) 加一个共享夹具头校验该层; ctest 全部 **7/7 通过**:

| 测试程序 | 用例 | 覆盖 |
|----------|------|------|
| `tests/test_manifest_reader.cpp` | 10 | 有效解析; 11 项不变量逐一拒绝; magic/版本检查 |
| `tests/test_mechanism_adapter.cpp` | 9 | 保持物种顺序; MW kg/mol → kg/kmol 换算; 全部 16 项 config 选项; `golaut` 与热泳护栏 |
| `tests/test_mom_context.cpp` | 9 | 从清单构造; `ComputeCell` 调用顺序; thermo/model 生命周期顺序; `diffusion_coefficient() > 0` |

`tests/test_manifest_fixture.hpp` 提供共享的 `make_test_manifest()` 辅助函数, 供三个程序构建一个有效的内存清单, 各测试再对其扰动以触发各拒绝路径。

### 16.5 构建说明

该层以 **MSVC 14.50** 通过 **Ninja** 生成器构建。与其余源文件一致, 施加 `/utf-8` 标志 (代码与清单含非 ASCII 字符)。**nlohmann/json v3.12.0** 单头依赖 **内置** 于 `mom/third_party/nlohmann/json.hpp` —— 无需网络拉取。所有新实现文件位于 `mom/src/MOM/`, 所有新测试位于 `mom/tests/`。
