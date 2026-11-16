# Loader-only V2 delta review

**Delta GREEN; actual attempt02 только после root loading-only control
и собственного resource admission.** Научные функции/domain/caps не
менялись. Reviewer Julia/HDF5/новые maps не запускал.

- Wrapper cb1d60a42c0be18761d5fa8e69846c1e5c8a81a08eb8b5c95ebc22483c31f144.
- False input49002f5cc19fbca991be36b7ba009b40eb6583a576b0eb0deada03bb8c0a1c73.
- Independently computed true-only input0057daf02b9c5fd80f0384350d6c18de35ac07795ed0369dc9045175459a8191.

Byte comparison с archived16512 подтвердил **единственную** source
замену после source-adapter include:

    Main.load_modules(sconf)
    → Base.invokelatest(() -> Main.load_modules(sconf))

Parsed input diff с archived7cfb — только wrapper_sha256. Refined ecf7/
MM9182/original1434/94c6, M3 qualification/review bindings и numeric/
OS/call/charge/cutoff settings unchanged. Старые true807183/input/log/
source сохранены в loader-v1-failed-before-model/SHA256SUMS.

## Почему это исправляет данный отказ

Actual attempt01 log показывает world-age warning **global binding**
Main.load_modules и MethodError нового метода, до source model/maps.
Closure вызывается в latest world; lookup Main.load_modules происходит
уже внутри неё, затем вызывается новый method. Простое
invokelatest(Main.load_modules,sconf) оставило бы lookup binding до
latest boundary, поэтому выбранная форма обоснована. Installed Julia
Base/reflection.jl:1342–1372 отдельно различает latest global resolution
и current-world argument evaluation; это согласуется с observed error.

Остальной путь проверен: Base.include — существующий Base method;
после load_modules выполняются includes Refined/MM, затем уже имеющийся
`Base.invokelatest(run_loaded,...)`. Все новые Main.source/plain/bytesout,
Refined pilot/MM correct и вложенные callback calls происходят внутри
этого latest-world входа, после завершения includes. Иных новых прямых
вызовов ещё не видимых bindings в этом startup path не найдено.

Root planned loading-only control должен проверить реальный include/
latest pattern в новой namespace без source HDF datasets/model/Φ.
Ранее13 tiny tests остаются evidence unchanged math/domain, но не
подменяют loader control. Новых raw/inner maps для повторения tests
не требуется. Actual failure1.82s/420276KiB,exit1 и более ранний missing-
input command failure сохраняются отдельно; оба0model/maps, не failed
charge/cutoff outcomes. Whole-file SHA preflight I/O не называется map.

Новый actual имеет те же prespecified512/1024/≤144/2GiB/1200s/1GiB
limits и обязан получить отдельный result review. Этот source delta
не даёт заранее numerical/scientific acceptance и не ослабляет guards.
