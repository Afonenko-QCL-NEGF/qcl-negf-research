# Pre-freeze review: finite-N / linear-M3 μ controller

**Actual HOLD сохраняется.** M3 API/module/qualification ещё не frozen,
M3_qualified=false/UNSET и runtime=false останавливают actual до HDF5.
Это draft review с ранними замечаниями, не final hash approval.
Reviewer Julia/Pkg/HDF5 не запускал, controller/source не редактировал.
Прочитанный current draft SHA2d50e3e326f586b995a9238c217dc602bab2378dcc0ab862cc4343235aaf385d,
202LOC; input wrapper_sha устарел относительно draft и должен быть
обновлён в финальном freeze. Source module1434/adapter94c6 неизменны.

## Ранние исправления до freeze

**F1, evidence preservation:** после MM.correct finite/eta/linear-residual
assertions стоят до density diagnostics/ρ write. При вычисленном,
но неквалифицированном M3 результате summary сохранит лишь текст ошибки,
без фактических η/condition/residual/number/Pauli и candidateρ. Нужно
записать candidate M3 witnesses с unaccepted label **до** этих gates;
только затем решать, допустимо ли предлагать density root value. При
исключении внутри MM.correct отсутствие result остаётся not_measured.

**F2, counters/failure stage:** maps увеличивается callback после
вычисленной Φ, а не до attempt. Это completed/callback evidence; strict
верхняя граница144 следует отдельно из≤8+1 calls×16 iterations. При
исключении до callback/возврата M.pilot текущий call ещё не добавлен в
calls.toml. До pilot нужно сохранить call-start/reserved-step record;
разделить started calls, completed callbacks и структурный upper bound,
не выдавать maps за точное число всех начатых maps. Retrying/extension
не требуются. Raw per-map snapshots уже сохраняются до causality gate,
history/calls пишутся до rejection unqualified inner — это правильно.

**F3, pending M3 qualification binding:** текущий boolean/hash gate
fail-closed, но при final freeze одной M3_qualified=true и SHA какого-либо
файла недостаточно для семантической квалификации. Нужен reviewed
qualification schema/status, exact M3 module/review/source/API bindings
и соответствующий positive result. Пока его schema неизвестна, final
approval невозможно; это не разрешение считать draft runnable.

## Что согласовано по текущему коду

- Source adapter/model guards выполняются при source μ. Новое context
  изменяет только μ и relative bare atom energies eps−Δμ; V/full H/K/wk/
  beta/Omega/thermal population shared. Это верное z=μ+iω convention,
  не перенастройка model/source μ. Ни Γ source, ни Poisson не читаются.
- Каждый inner отΣ0,16steps,α1. Status iteration_budget и last3 raw
  REL≤1e−9 **AND** ABS≤1e−12, finite, нужны до предложения scalar density.
  Norm — unweighted Frobenius по n,k,a,b в nativeΣ units; wk-weighted
  numerator/scale отдельно diagnostic, без β/charge renormalization.
- Native target .045/sourceE0 .1eV: bracket±.001native означает±.1meV.
  Sign bracket/secant с midpoint fallback не требует monotonicity;
  нет expansion. Strict interior guard ловит exhausted Float spacing.
  Все8 N128 calls включают два endpoints; максимум1 N256 call на
  **той же** selected μ, тоже16steps отΣ0, ограничивает maps144.
- Финальные gates проверяют corrected **number** в±1e−8target и
  N128/N256 **number difference**≤1e−8target. Это не matrix-density
  cutoff bound. Rawρ/one-shot/linearρ и Pauli сохраняются без clipping;
  selected number может отличаться от raw finite-N number по новой
  M3 closure. Tail/infinite-N/SCBA certificate не следует.
- Stage использует frozen source model reconstruction, затем small
  Mats arrays; сохранённый native full-cache planning budget source
  не даёт права native solve в OS2GiB. Source adapter old workload gate
  обходится явно, воспроизведены source/env guards и preserved source
  checktime; новый wrapper 2GiB/1200s/J1/BLAS1/CPU≤200% задаёт собственный
  scope. OS enforcement не заменяется heap hint.
- Global output512MiB/file32MiB,1MiB final reserve, отдельная directory
  каждого μ/N, immutable per-step S/G/ρ и latestΦ. Peak/стоимость ещё
  not_measured; callbacks сохраняют raw current maps, не смешанные steps.

Неблокирующее улучшение: endpoint b вычисляется до проверки, что endpoint
a уже достиг charge target; это может расходовать ненужный call и
потерять найденный endpoint при failure b. Можно проверить a раньше,
не расширяя budget. Matrix-density norm witness между selected128/256
также полезен, если пользователь хочет вывод оρ, а не только scalar N;
нынешний number-only gate надо так и называть.

После F1/F2 и frozen M3 API/F3 нужна новая source/hash delta review и
root source-free controls, затем отдельное actual resource decision.
Вопрос о корректности самой M3 asymptotic/linear closure решается её
отдельными доказательством и qualification; root/controller не принимает
её задним числом по одному successful charge root.

## Исправленный μ128/256 draft, затем решение UNRUN

Reviewed e074a1f726500cbeff3fb92a59dc085c96f10740ff16c2df5333416315cd3dc8
и последующее exact8source-key tightening05c9a86c90407a69033237d494ae05cb2e2bb577d4ddcc43c84e5c3ad56f3c01
(input80499bc8…false). F1 закрыт для returned MM.correct results:
candidateρ/S0/B0 и full diagnostics сохраняются unaccepted до caller
finite/eta/residual gates; positive scalar-proposal qualification отдельным
файлом. Исключение внутри MM.correct остаётся explicitly not_measured
после M3-request, не подставляется нулевой результат.

F2 закрыт: request/call status started_reserved и16-step reservation
пишутся до pilot; completed callbacks и reserved/structural144 разделены.
F3 implementation исправлена: actual postprocess complete/zero maps/μfits,
BLAS initial/final1, q/qi text/source equality, module/adapter/source μ/
config/source-file digests и положительные finite N128/N256 linear-solve
witnesses проверяются. Источник содержит exactly8 required files;
пустая qualification table больше не проходит vacuously. Rebuilt model
binding сравнивается exact с qualified source-binding. Explicit BLAS setter
и getter до работы/в re-entry устраняют надежду только на environment.
Финальный M3 review/input freeze всё равно обязателен; actual не принят.

Root затем решил **не запускать μ128/256**: actual M3-source corrected
N256−N128=−9.172036394e−9, т.е.2.03823e−7target,≈20.38× planned1e−8
cutoff. Это не доказательство такого mismatch на ещё не найденной μ,
но уже не основание предполагать данный cutoff достаточным. Старые
draft/source/input сохраняются UNRUN; это decision по measured evidence,
не failed μ experiment.

Новый higher-cutoff candidate outer512/witness1024 — отдельная identity,
пока не reviewed здесь. Нужны freeze и independent checks module-copy
body equivalence, N-generic counters/gates/output envelope, а также
M3 result qualification. Current-HOLD не снимается историческим closure
F1/F2/F3 и lower-N postprocess success.
