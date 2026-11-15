# Независимая классификация warm attempt01 obstruction

**Warm attempt01 остановлен на strict spectral-weight domain guard до первой Φ.** Это numerical/domain obstruction выбранного initial retarded evaluation, не физическое доказательство отсутствия равновесного решения, не SCBA convergence/failure и не новый tolerances pass. Guard сработал; warm16 green был условен прохождению initial gate, который здесь не достигнут.

Primary `equilibrium-gamma-warm/attempt01`: summaryexecutionstopped,π_calls1/precheck_actual_calls1, всеnative events/start/completion bounds0,small_checks pass. Logexit1,39.97s processwall,26.0509s harness,peak3288032KiB,swap0. Нет OOM/timeout. Model/Gamma-first/seed-manifest/warm-initialization сохранены; нет map001, actual-one-map, Φ checkpoints либоnative result. Рецензент читал только TOML/log/code, не большие arrays и не выполнял Julia/Pkg/maps.

`warm-initialization.toml` подтверждает manifest18439792/seed4142e06d/sourcea690 и completed actual rebuilt-problem fingerprints: grids/profiles/kernels/qK/basis/h соответствуют cold target. ImportedΓcone min+1.25852858e−27, blockF3.62865e−27, Hermitian defect0 при backwardallowance7.1054e−14. Таким образом отказ не является незамеченной подменой seed/frame и не произошёл в initialization/source guard. Малое положительное computedΓ не является строгим PD certificate, но evidence явного eigen-clipping нет.

Точное место: `evaluate` строитretarded_state, A, пропускает A через declared backwardcone, получаетq из owning `_seed_number_weights_production`, затем `fit_mu` падает на `all(isfinite,q)&&all(>=(0),q)`. Zero разрешён: сообщение «nonpositive/nonfinite» шире фактического условия. AcceptedCore helper сам отвергаетnonfinite доreturn. При неизменном source/control flow эта остановка устанавливает наличие **отрицательного конечного q_e**. Значенияminq/count/location/negative mass/error scale этими artifacts ещё не измерены. Нельзя объявить их «только harmless roundoff» без следующего witness.

Nativeq =gs*wE/(2π)Σk wk tr[i(GR−GA)], положительные weights; отрицательныйq означает отрицательный вычисленный summed spectral trace на некоторых energy nodes. Acone допускает absolute floating backwarderror и поэтому логически может пройти при отрицательномq, тогда как μ monotonicity/capacity fit требует rawq≥0. Ослаблять последнее условие или подменятьq=max(q,0), новуюη или Γfloor нельзя по факту этого отказа.

Следующий ограниченный diagnostic должен различить минимум:

- retarded-field inconsistency: Γinput против i(ΣR−ΣR†), anti-Hermitian contamination HPI иcavity/totalΓR, относительно actual tinyinput broadening;
- inverse/subtraction arithmetic: native GR/Acancel=i(GR−GA) против stable resolvent identity и independently constructed central45 solve на **том же** declared operator;
- существующий model/input noncausal part или boundary issue, если ошибка уже присутствует в precise declared D и не исчезает при increased arithmetic precision.

Для BigFloat witness нужно явно различить high-precision solve **float64D как exact input** и пересчёт ideal HPI/Γ model: первое изолирует inverse/subtraction, но наследует округлённый/possibly non-Hermitian ΣR. Нельзя назвать его всей physical моделью. Selectedbadnodes должны иметьq/minblock scales, condition/backward residual, Γ/HPI/cavity matrices и units; aggregateq используетallk, поэтому единичный selectedk не воспроизводит wholeq автоматически.

Π не дошла доμfit success, FD Green, scatteringΦ или initialKMS comparison. Currentμ/rΓ/AA trajectory иnative audits **not_measured**, oldsource μ иoldresidual не переносить. Правильное решение — сохранить этот отрицательный execution/domain outcome и разобрать отдельный frozen diagnostic; не продолжать warm16 в обход guard. Static diagnostic review ожидает exactfreeze.
