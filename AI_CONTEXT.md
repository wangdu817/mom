# AI Collaboration & Architecture Rules for libDictionary

## 1. Technical Stack, Standards & Style
- **Language Level:** Pure Modern C++ (strictly targeting C++20 and C++23 standards).
- **Function Declarations:** Prefer trailing return types (`auto func(...) -> Type`) for complex or highly templated functions to improve readability.
- **Attributes & Keywords:** Enforce code correctness and intent at compile time:
  - Mark all non-void, non-side-effect functions with `[[nodiscard]]`.
  - Apply `noexcept` to all functions with strong exception-safety guarantees (especially move operations, getters, and core mathematical utilities).
  - Use `constexpr` aggressively for arithmetic and configuration evaluators.
- **Dependency Philosophy:** Maintain a minimalist dependency footprint. Core utilities must rely *exclusively* on the C++ Standard Library. Use Boost C++ if really convenient. Heavy ecosystems must remain entirely optional and isolated.

## 2. Library Distribution & Compilation Strategy
- **Binary/Header Split Architecture:** The library is designed as a hybrid pre-compiled structure to protect proprietary source code (`.cpp` implementation files) while maximizing user usability:
  - **Public Space (`include/`):** Contains interface definitions (`.hpp`) and template implementation bodies (`.tpp`).
  - **Private Space (`src/`):** Contains concrete implementations, adapter bindings, and explicit template instantiations.
- **Template Hiding:** Use `extern template` declarations in the public headers to prevent the compiler from implicitly instantiations templates in user translation units, forcing the linker to bind against the pre-compiled static/shared library binary instead. Ensure AI suggestions do not break this boundary.

## 4. Performance & Hardware Optimization
- **Data Layout:** Prioritize cache-locality. Prefer contiguous memory structures (like `std::vector`) over linked or fragmented structures.
- **Memory Overhead:** Eliminate redundant memory allocation loops. Use `std::string_view` for read-only string tracking, pass heavy objects by `const&`, and leverage move semantics (`std::move`) natively.

## 5. AI Interaction & Output Preferences
- **Conciseness:** Provide zero-fluff, highly direct C++ code modifications. Omit lengthy introductory or concluding conversational prose.
- **Rationale Over Guesswork:** When suggesting code optimizations, explicitly back the decision with a low-level hardware rationale (e.g., SIMD vectorization friendliness, cache line bouncing prevention, or stripping down stack-unwinding plumbing via `noexcept`).
- **Build Awareness:** Understand that this project uses CMake 4+ and maps targets cleanly via modern target-based dependency properties (`target_link_libraries`, `target_include_directories`). Do not suggest legacy global variables like `include_directories()`.