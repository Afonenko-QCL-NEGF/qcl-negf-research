# Независимая проверка final root pair

Решение: **APPROVE фактически выполненный технический integration smoke** пары Core `db26a0ef49673d8246ae336cbe73bafe6d7482ac` / Runner `53b3a280109228a980df7ddea59e4f1b75b836be` в path-developed Pkg environment. Научной приёмки, приёмки Runner checkpoint/restart на новой паре или standalone pin на db26 этот опыт не даёт. Перед окончательным handoff нужно исправить устаревший reproducer в `environment/final-integration/HANDOFF.md`, как указано ниже.

## Исходники и область решения

- Проверены реальные HEAD/tree и пустой `git status --porcelain` обоих detached worktree: `core-final-integration` tree `073abe5070b6ca60822f5a9b1b139c614ec3e16f`; `runner-final-integration` tree `3d6ec2b9372829ebfd65b6e98265e8ad590a99b5`.
- Core diff от принятого научного `e171caa60cce676c6c4238dd671e7c527f2bd977` имеет SHA256 `14d3be2d95786d0ad5c5e458e99e24dd8ebcc8ee183a812efdc8f7e92dd72a8a`; `src/test` — `8b47d71a2beedee65dfd7cb1620ebc8baa704870db3d9b1fb16d98d6d7cf2d7e`; четыре src — `00e2e7fad176b9411794b8c756f39406b850281b870ec801d83d01658593f751`. Это точные ранее проверенные IFR publish/code/source patches, без дополнительной дельты. Runner tree идентичен ранее проверенному `757d8cf5ebeaba013ff4b30235327f9dc75ae404`.
- `core-accepted-science` остаётся чистым на e171/tree `742ed5a6fe82f1e5c63fc1a2290ed77f89ca6da8`. Результаты кампании на нём не переименовываются в результаты db26.
- На момент проверки root HEAD `b515c42ef14c1f61151626280da0bcf100d116be`, рабочие component checkouts Core f50/Runner18f и изменённые gitlinks отличаются от данной подготовленной пары. Это свидетельство для последующего явного обновления root gitlinks, а не доказательство уже опубликованного состава root.

## Pkg и runtime evidence

`environment/final-integration/prepare-root.log`: штатный `tools/prepare_workspace.jl ... --runtime-only`, native `Pkg.develop` обеих компонент и `Pkg.instantiate`, offline/precompile-auto=0; exit0, 7.83 s, maximum RSS515124 KiB, swap0. Эффективные Project/Manifest содержат относительные path sources именно на проверенные worktree; Manifest Julia1.13.3.

`run-check.sh` сначала штатно добавляет **прямые зависимости checker** Unitful1.28.0 и HDF5 0.17.4. Они уже присутствовали транзитивно после prepare; `check-root.log` сообщает `Manifest No packages added to or removed`. При независимом сравнении с `environment/ifr/Manifest.toml` все общие dependency entries, кроме намеренно отличающихся Core/Runner, равны; добавлены YAML и StringEncodings для Runner. Это не универсальное утверждение о неизменности всех исторических окружений.

`check-root.log` подтверждает оба фактических import paths, исполнение assertions точных HEAD и `FinalPairChecksOK`; exit0, 45.86 s (включая Pkg.add), RSS903096 KiB, swap0. Native command: Julia1.13.3, threads2, BLAS1, bounds checks, offline/precompile-auto=0. Heap hint сам по себе не является hard cap; фактический RSS находится ниже обсуждаемого resource window.

## Что именно проверил checker

`check_pair.jl:1–41` импортирует Core, Runner, Unitful и HDF5; проверяет `pkgdir` и HEAD. Через Core public API строит reference parameters с tutorial grid NE49/Nk5/Nb3/Nz48, acoustic-only scattering, `AlgorithmOptions(kernel_build=:direct)`; остальные алгоритмы остаются default `hilbert=:fft`, `elastic_retarded=:hilbert`. Выполняет Core `solve_scba_production` при фиксированном U=0 и max_scba=2; Poisson loop не вызывается.

Проверены две строки history и finiteness полей r_D/r_A/r_K/r_Σ/r_λ/λ, затем точное чтение real/imag GR/Gl после прямой записи `HDF5.h5open`. Ни Runner execution/serializer/loader, ни новая IFR contraction здесь не вызываются (IFR=false). Следовательно, это совместный импорт пакетов плюс Core solve и direct HDF5 roundtrip; не Runner checkpoint/restart roundtrip и не повтор IFR benchmark.

Явный результат: `FixedUStatus=max_iterations converged=false iterations=2`; существенные residuals второй итерации rK≈0.958, rΣ≈0.991. Finiteness/serialization и малый rD не дают научной сходимости, grid adequacy или экспериментальной валидации. Новый runtime reviewer не запускал; HDF5/массивы не читал.

Отдельная ранее принятая standalone Runner проверка 21 loader/checkpoint assertions относилась к его immutable source pin **Core e171**. Runner Project/docs до сих пор содержат этот pin. В данной integration env он штатно перекрыт явным Pkg.develop Core db26; эти два маршрута имеют разные source identities и их свидетельства не объединяются в новую standalone checkpoint приёмку db26.

## Reproducer gap и требуемая правка handoff

Проверенный `run-check.sh` содержит необходимые две `Pkg.add(PackageSpec(...))`, но `HANDOFF.md` в проверенной версии SHA `f661502dbefdcfcb15ffc5e3ee4bffbf422f24c584fa185f4787f4cafe31c0bb` по-прежнему говорит «не выполнялись» и предлагает после runtime-only prepare немедленно запускать checker. Runtime-only Project с двумя пакетами не делает их транзитивные Unitful/HDF5 доступными как direct imports main project.

В окончательном handoff нужно отразить выполненные логи/refs и добавить этот Pkg.add (либо сослаться на фактически выполненный run-check.sh после prepare). Checker намеренно отвергает существующий output; сохранённый `fixed-u-smoke.h5` не удалять ради повторения, использовать новый отсутствующий output path/копию script с новым назначением. Также сохранить native resource envelope и ENV фактической команды. Эта документационная поправка не требует нового Julia runtime; источники и успешный actual route уже проверены. Reviewed source/repo/index не менялись.

## SHA256 первичных малых артефактов

| Файл в environment/final-integration | SHA256 |
|---|---|
| Project.toml | `2884e04c31ab2968f06aa75c8f070116816a1f2d15a47e766eab65a70fb87151` |
| Manifest.toml | `27bd15346d6f9de90427f643ecbba14ac7785f1074a270b8f2519e894b90d7a8` |
| check_pair.jl | `1b7dfb5b930b6d0fd1a74cb1989919180bb067f26e7204e412f848972ac45969` |
| run-check.sh | `9d857d64bd6aa46b839d2a671bb0a4849c7418054d6ca82cbf7bc59e1a8b92da` |
| prepare-root.log | `ee3e42b2c2f6302f843fc0945c1aa53a28dd11964665ad315a306730d8bace5c` |
| check-root.log | `8d7be398a18e664a39ab43fc3062cdd3ddd8812ecd8306253342e130021f4efd` |

Существующие owning Core/Runner reviews и `exact-ifr-benchmark-result.md` остаются основаниями для соответствующих отдельных контрактов; данный smoke связывает конкретную будущую root пару, не расширяя их научных claims.

## Handoff closure

Root исправил `HANDOFF.md`; независимо прочитанная новая версия имеет SHA256 `f0b333e2d826d3684db700a93487622197ad279a71d4c0ab52fdc45dbbb94ded`. Теперь указаны actual45.86s/exit0, direct Unitful1.28/HDF50.17.4 Pkg.add, новый отсутствующий output path и точная narrow scope. Документационное замечание закрыто; исходники и выполненный checker не менялись, новый runtime не требуется. Reviewer handoff не редактировал.
