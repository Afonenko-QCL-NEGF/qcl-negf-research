# Маршрут исследователя

Research хранит цель, модельные inputs, варианты и источники. Runner разрешает их в effective
configuration и frozen plan; Core владеет уравнениями и численными критериями. AiiDA связывает
execution/attempt с исполнением, Results читает сохранённые данные, Portal представляет те же
публичные операции. Подготовка plan, чтение результата, явный derive и новый solve — разные
действия.

## От вопроса к используемым значениям

Сначала выберите вопрос, наблюдаемую и различающее свидетельство. Укажите предпосылки, граничные
условия, единицы и критерии применимости. Отделяйте проверку реализации от достаточности сетки и
сравнения с экспериментом. Для equilibrium route используйте [S01 brief](equilibrium-s01.md); он
сохраняет неизвестные origins и открытые критерии, не объявляет candidate принятым.

Study задаёт `configuration.sources` в порядке merge, затем собственные `configuration.overrides`,
варианты и axes. Их смысл разрешает один Runner; не дублируйте material resolver или формулы в форме
UI. Перед запуском проверьте resolved preview, point axes, protocol/material/model identities и
variant diff. CLI `configure` и его bounded JSON API описаны в
[owning Runner configurator](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl/blob/codex/project-simplification-2026-10-27/docs/src/user/configurator.md).
Используйте API совместимой установленной версии, а не самостоятельно составленный
неверсионированный JSON packet.

| Вид происхождения    | Что устанавливает                                        | Чего недостаточно                                                    |
| -------------------- | -------------------------------------------------------- | -------------------------------------------------------------------- |
| Formula/model origin | Источник уравнения, approximation и области применимости | DOI о структуре не обосновывает все material constants               |
| Value origin         | Обоснование численного значения с единицами и условиями  | Последний YAML assignment сам по себе не литературный источник       |
| Assigning origin     | Какой source/override последним назначил leaf            | Отсутствующая leaf history не доказывает происхождение typed default |
| Used snapshot        | Effective значения конкретного выбранного point          | Текущий checkout и начальный preset не заменяют frozen inputs run    |

Новый override создаёт явно выбранный вариант постановки. Не уменьшайте scattering, сетку, базис,
окно или допуски ради скорости/сходимости. Unknown source допускается как явно обозначенная
модельная предпосылка. Из `x_Al` нельзя выводить новые параметры материала без owning model rule.

## Проверки до исполнения

[Структурная проверка](../README.md#структурная-проверка-inputs) не вызывает Core. Затем совместимый
Runner может разрешить выбранный вход:

```sh
mkdir -p plans
qcl-negf plan config/equilibrium.yaml > plans/equilibrium.json
```

Проверьте, что это один independent point 0 mV/period, 70 K и baseline variant. Это **план будущего
S01**, а не выполненный научный опыт. Frozen plan содержит resolved inputs и fingerprint; исходные
source files и model passport объясняют выбор. Наблюдаемые значения при adaptive domain могут
требовать отдельного actual state context, который не подменяет исходный model snapshot.

Resource estimate — прогноз с указанными unknown/confidence; он не является выданной памятью.
Platform/Slurm устанавливают actual resources и admission. Перед исполнением нужны свежие данные о
host reserves, storage, runtime source и состоянии worker. Общих campaign quotas по
CPU/RAM/времени/attempts этот маршрут не вводит. Конкретный job использует согласованный resource
envelope и stop rules для нарушений инвариантов; нехватка capacity не разрешает менять физическую
модель.

## Новый запуск выбирается отдельно

Инженерная приёмка системы использует synthetic jobs. После неё отдельное научное решение выбирает
S01 и критерии. Только на этом этапе команда Runner может выполнить frozen plan:

```sh
qcl-negf run-plan plans/equilibrium.json results/equilibrium
```

Это дорогой stationary solve; при работе с кластером тот же plan передают через публичный AiiDA
submission с зарегистрированным Code и actual Slurm resources. YAML научной постановки не выделяет
ресурсы. Continuation принадлежит Runner внутри execution; scientific failure не становится accepted
seed. Независимые branches могут продолжаться по своей политике. Exact аварийный resume сохраняет
постановку и не означает автоматический retry научной неудачи.

После mutation timeout сначала сверяют прежний request/attempt и receipt. Не отправляйте ту же
операцию вслепую с новой identity. Отключение клиента само по себе не отменяет работу.

## Чтение, сравнение и отчёт

Выберите точные run/execution/attempt, source/model/plan и completeness. Results читает capabilities
сохранённого формата; missing capability должна иметь причину. GET, открытие help и построение
графика не запускают solve или повторную scientific acceptance. Derive — отдельная явно выбранная
операция с указанными prerequisites; недостающие данные не заменяют нулями.

Для сравнения вариантов проверьте сопоставимость структуры, bias/temperature, scattering,
grids/basis/windows, normalizations и criteria. Большой payload хранится один раз в научном storage;
AiiDA связывает малые locators/identity/status и provenance. Backup и метаданные не являются лишними
научными копиями.

Отчёт ссылается на точный state/artifact и отдельно отвечает на пять вопросов: выполнился ли
процесс, сошлись ли уравнения, пройдены ли физические проверки, достаточна ли дискретизация, есть ли
экспериментальная validation. `RUNNING`, exit0 и совпадение двух несошедшихся candidates не дают
положительных ответов на остальные вопросы. `not_measured` остаётся неизвестным. Для near-zero
current нужны обоснованные абсолютные единицы/ошибка; относительное деление на ноль не является
критерием.

| Утверждение                    | Свидетельство                                     | Ограничение                                              | Решение                |
| ------------------------------ | ------------------------------------------------- | -------------------------------------------------------- | ---------------------- |
| Постановка S01 доступна        | Versioned study и source/model passport           | Scientific dataset ещё не предъявлен                     | сохранить              |
| Plan можно проверить без solve | Public Runner `plan`, structural validator        | Matching runtime и semantic receipt проверяются отдельно | дополнительно измерить |
| Candidate научно принят        | Требуется отдельный опыт и перечисленные evidence | Инженерная приёмка этого не доказывает                   | данных недостаточно    |
