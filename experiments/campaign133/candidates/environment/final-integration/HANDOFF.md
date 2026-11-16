# Фактически проверенный прямой запуск root-пары

Core `db26a0ef49673d8246ae336cbe73bafe6d7482ac` и Runner `53b3a280109228a980df7ddea59e4f1b75b836be` импортированы из отдельных чистых detached worktrees. Штатный offline Julia Pkg подготовил окружение; для checker явно добавлены Unitful1.28.0 и HDF50.17.4, уже присутствовавшие транзитивно в Manifest. Глобальное окружение не менялось.

`run-check.sh` сохраняет исполненную команду. Pkg preparation:7.83s/515124KiB. Добавление direct checker dependencies и проверка:45.86s/903096KiB/swap0/exit0. Checker подтвердил actual paths/HEAD, две SCBA итерации при фиксированном U=0 и HDF5 roundtrip массивов. Статус max_iterations/converged=false. Это технический smoke, не научная приёмка, не Runner execution/checkpoint contract и не упражнение нового IFR operator(acoustic-only model). Отдельные прежние Runner checkpoint tests относятся к Coree171.

Для нового воспроизведения из корня superproject после подготовки этих exact worktrees, с новым output path:

```sh
task_env="$PWD/.build/campaign133/environment/final-integration"
task_core="$PWD/.build/campaign133/core-final-integration"
task_runner="$PWD/.build/campaign133/runner-final-integration"
JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0 julia --startup-file=no "$task_runner/tools/prepare_workspace.jl" "$task_core" "$task_runner" "$task_env" --runtime-only
JULIA_PKG_OFFLINE=true JULIA_PKG_PRECOMPILE_AUTO=0 julia --startup-file=no --project="$task_env" -e 'using Pkg; Pkg.add([Pkg.PackageSpec(name="Unitful",version="1.28.0"),Pkg.PackageSpec(name="HDF5",version="0.17.4")])'
OPENBLAS_NUM_THREADS=1 julia --startup-file=no --threads=2 --check-bounds=yes --project="$task_env" "$task_env/check_pair.jl" "$task_core" "$task_runner" "$task_env/new-fixed-u-smoke.h5"
```

Root ограничивал процесс3GiB/high2.5GiB/swap0/CPU200%,180s,file64MiB. Перед новым запуском проверить ресурсы, exact refs/clean state и отсутствие нового output. Initial prepare-only handoff сохранён в HANDOFF.pre-run.md; source/log/Project/Manifest/checker сохранены без изменения результатов. Научные опыты кампании остаются на frozen Coree171.
