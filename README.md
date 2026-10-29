# QCL-NEGF Research

Версионируемые научные inputs и литературные источники для QCL-NEGF. Это репозиторий постановок;
готовых результатов расчёта здесь нет. Начните с
[русского руководства исследователя](docs/researcher-guide.md).

Уравнения и вычисления принадлежат [Core](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl),
разрешение конфигурации и CLI — [Runner](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl),
provenance и Slurm workflow — [AiiDA](https://github.com/Afonenko-QCL-NEGF/qcl-negf-aiida). Root
gitlinks выбирают совместимую композицию компонентов.

| Вход                      | Научная цель                                                                                           | Выполнение                                                                          |
| ------------------------- | ------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------- |
| `config/operators.yaml`   | Семь алгебраических, метрических, континуальных, cavity, LO equilibrium и spectral resolution проверок | Диагностические вычисления                                                          |
| `config/spectra.yaml`     | Спектры структуры и сравнение базиса для литературной структуры и синтетической ямы                    | Диагностические вычисления                                                          |
| `config/equilibrium.yaml` | Один независимый S01 candidate при 0 mV/period и 70 K                                                  | Дорогой stationary solve; [постановка и открытые критерии](docs/equilibrium-s01.md) |
| `config/convergence.yaml` | Независимый candidate при 48 mV и девять контролируемых численных сравнений при 70 K                   | Дорогие stationary solves                                                           |

Все включённые studies достижимы из этих четырёх входов. Convergence bundle содержит 30 721 energy
nodes для reference candidate и до 46 081 для сравнения; это не smoke test ноутбука.

## Планирование без solve

В подготовленной совместимой установке `qcl-negf`, из этого репозитория:

```sh
mkdir -p plans
qcl-negf plan config/operators.yaml > plans/operators.json
```

`plan` разрешает inputs и создаёт frozen JSON plan с fingerprint. Он не выполняет transport solve.
Сохраните выбранный plan вместе с будущим run. Проверка planning не доказывает наличие ресурсов,
сходимость или физическую корректность результата. Новый научный запуск выбирается отдельно после
инженерной приёмки; команда запуска описана в руководстве, её не выполняют для проверки
документации.

## Структурная проверка inputs

Используйте Python 3.14 и зависимости, заранее подготовленные для выбранной Root-композиции. Из Root
checkout:

```sh
uv run --locked --no-sync python components/qcl-negf-research/tools/validate.py
uv run --locked --no-sync pytest components/qcl-negf-research/tests
```

`--no-sync` требует уже подготовленное окружение. Его установку выполняют отдельным шагом по
текущему Root runbook. Validator проверяет граф studies, merged variants по общим schemas, duplicate
YAML keys, missing/escaping paths, недостижимые fragments и контрольные суммы литературных данных.
Он не запускает solver. Семантическая проверка `qcl-negf plan` — отдельная проверка Runner.

## Интерпретация и изменение постановки

- [Руководство: цель, inputs, used values, run и чтение](docs/researcher-guide.md)
- [S01: thermal equilibrium candidate и открытые научные gates](docs/equilibrium-s01.md)
- [Model assumptions и контролируемые сравнения](docs/scientific-inputs.md)
- [Литература, provenance и пределы сравнения](literature/README.md)
- [Правила вклада и проверки](CONTRIBUTING.md)

Код и документация используют [MIT license](LICENSE). Метаданные публикации и оцифрованные данные
сохраняют атрибуцию в
[literature/reference-2019/source.yaml](literature/reference-2019/source.yaml); статья не
распространяется в этом репозитории.
