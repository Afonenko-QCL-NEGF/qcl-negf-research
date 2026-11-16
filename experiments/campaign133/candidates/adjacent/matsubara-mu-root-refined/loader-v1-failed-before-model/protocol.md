# Передача refined μ-root wrapper: WRITE-ONLY/HOLD

Автор не запускал Julia/Pkg, source HDF5 или μ experiments. Новые файлы
ограничены этой новой папкой; original module1434/sourceadapter94c6 неизменны.
Незапущенный draft2d50 и его SHA сохранены в `../matsubara-mu-root/unexecuted-draft-v1/` перед
исправлениями BLAS/F1/F2/F3 по review. Новый refined draft содержит231LOC; предыдущий05c9/input8049 UNRUN сохранён в ../matsubara-mu-root/unexecuted-refinement-predecessor/.
Дизайн и input задают отдельную linear-M3 fixedcharge identity. Текущий
input=false запрещает runtime; M3_qualified=true теперь означает **только
принятый метод postprocess** по result review72598398…, не charge/cutoff/
infinity pass. M3 N128/256 actual cb2dcc… завершён, все bindings заполнены.
Предыдущий input5c86 сохранён как input-prebinding-unexecuted.toml. Заполненный
input7cfb6195… ожидает отдельного binding review; Root вправе не запускать его.

M3 module/API, tiny и actual postprocess уже имеют frozen свидетельства.
Новый refined domain/runtime требует собственного static/binding review,
Root tiny и отдельного решения об actual. Никакого
автоматического fallback, tail-order изменения или смены source identity.

Только после static review Root может выполнить source-free tiny (без H5/maps):

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 julia --startup-file=no --project=.build/campaign133/environment/accepted -e 'include(".build/campaign133/adjacent/matsubara-mu-root-refined/run_mu_root.jl"); @assert BoundedMuRootRefined.sourcefree_controls().pass'
```

Wrapper явно вызывает BLAS.set_num_threads(1) при старте и перед preflight,
source build/inner pilot; getter1 проверяется и сохраняется. ENV переменные
сами по себе не ограничили установленный OpenMP OpenBLAS (Root witness128).
Это проверяет REL+ABS last3, Inf/NaN rejection, context μ/pole shift и secant
fallback. Tiny не квалифицирует физику или runtime. Новая scalar root схема
получает sign bracket и не предполагает монотонность плотности.

Actual возможен только в новом ROOT-owned output и explicit authorized copy
замороженного input; меняется только runtime_authorized. Root hard-envelope:
2GiB/swap0,1200s на весь процесс, Julia1/BLAS1, CPU≤200%, FS32MiB,
data1GiB. Не более8 N512 calls+одного N1024, каждый16 maps изΣ0/α1.
Wrapper не вызывает original adapter gate, чей input envelope относился к
старому опыту; source guards воспроизведены перед original load_modules/source
при sourceμ. Внутренние source build/checktime остаются неизменны.

```text
julia --startup-file=no --project=.build/campaign133/environment/accepted run_mu_root.jl AUTHORIZED_INPUT NEW_OUTPUT
```

Сохранение: каждый μ/N имеет отдельную папку, каждый map rawΣ/G/finiteNρ,
latestΦ и raw scalar norms; wk-weighted raw — только диагностический.
Terminal: raw/one-shot/linearρ, числа, разность closure, t4 bounds, linear
residual/condition/η и final Pauli. Транзиентные отрицательные eigenvalues
rawρ сохраняются без clipping; никаких Pauli-pass всех итераций не объявляется.
M3 candidateρ/diagnostics пишутся с candidate_accepted=false **до** finite/η/
linear-residual guards; positive qualification — отдельный малый TOML.
При исключении внутри correct M3-request остаётся not_measured. Перед pilot
call-start/reserved16 записаны; maps_completed_callbacks — только завершённые
callback, не число всех начатых вычислений. Structural upper bound144 отдельно.
M3 gate parse требует completed_postprocess_only/zero inner maps/μfits,
actual BLAS1, module+adapter SHA, source/input tables, два finite qualified
linear solves и equality qualification model binding с source context.
1MiB output reserve оставлен итоговому fail summary; hard-kill может оставить
tmp и последнюю сохранённую карту, не получает скрытого success.

`qualified_fixedcharge_M3_variant` означает только finite-N raw qualification,
charge/cutoff checks этой closure. `native_pass=false`,
`scientific_acceptance=false`, `infinite_N_certified=false` остаются явно.
Residual исходного finite-N map и M3 linear residual — разные поля; ни один
не является оценкой полной infinity-tail ошибки.

Refined module ecf7b1cc…: only name+allowedN+comment differ from1434; numerical bodies unchanged. Original source(c_source) still1434 frame, MM.applyK original identical. Outer charge1e−9; independent witness charge/cutoff1e−8. Blas explicit setter/assert and Root OMP1. Actual M3 N128/256 report qualification certifies this postprocessor execution only, never N512/N1024 accuracy.
