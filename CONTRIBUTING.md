# Изменение научных inputs

Создавайте study в `config/studies/` с явными configuration sources, единицами, научной целью и
контролируемыми вариантами. Включайте его в подходящий top-level plan. Все fragments должны
оставаться достижимыми; общие model parameters не копируйте между studies.

Для исправления текста/ссылок проверяйте diff, изменённые ссылки, anchors и форматирование. Julia,
установка зависимостей, docs build и transport solve для docs-only не требуются.

При изменении study inputs или validation code выполните owning schema checks и tests из
[README](README.md#структурная-проверка-inputs) в заранее подготовленном matching окружении.
Семантику каждого изменённого entry point проверяет `qcl-negf plan` совместимой версии Runner. Это
отдельная Julia-проверка без transport solve; подготовка окружения и научный запуск — отдельные
действия. Отсутствующий runtime ограничивает свидетельство, не блокирует независимое исправление
текста.

Изменение assumptions, layer dimensions, scattering или interpretations критериев требует научного
обоснования. Ссылайтесь на публикации по названию и DOI; сохраняйте единицы, digitization
uncertainty и условия. Не меняйте label непроверенного reference candidate, чтобы выдать его за
accepted result. Свежие ресурсы, storage reserves и actual admission определяют возможность
исполнения; общих campaign quotas по CPU/time/RAM/attempts этот документ не вводит.

Результаты храните в научной data/provenance системе. Machine-specific paths, большие массивы,
credentials и execution logs не коммитьте сюда. Root gitlinks определяют состав; интеграционные
проверки выполняются из [qcl-negf](https://github.com/Afonenko-QCL-NEGF/qcl-negf). Untrusted
contributions проверяются на изолированных disposable runners.

Рабочий маршрут и различие formula/value/assigning/used origins изложены в
[руководстве исследователя](docs/researcher-guide.md). Новое исследование проверяет свои критерии
отдельно от инженерной synthetic-приёмки системы.
