# Stable-source analytic Jv: независимая проверка результата

**Решение: bounded diagnostic принят как локальное свидетельство derivative consistency и ограниченного Krylov response.** Analytic covariance Jv согласуется с central FD второго порядка в направлении v=Γ. Arnoldi3 не даёт надёжной моды около1 и предсказывает лишь1.59% уменьшения restricted linear residual; причин плато или успешного Newton шага не установлено. Γ не обновлялась, SCBA/native acceptance остаётсяfalse.

Actual `adjacent/stable-source-jv/attempt01`, summarySHA da9ff0bda035851ef1dd53ef61cf73cee405ab5088d5b0b4982f277a84bc8391, runtimeinput1d3ed9d80db7acf50337afba75a0057d7b2f117883d9c96e0cd8aaf6c64849fa, adapter6beaeafe…/module5a56b333…/serializer6157cdb7…/stream8df…/parent50bc…/spectraleeca…. Все архивированные source hashes перепроверены по actual identity; runtime input отличается от frozen толькоauthorized bool. Coree171/tree742, Julia1.13.3/J4/BLAS1, fixed hH/HermPI/V0/finite4/model700ae…/Γlast64cab810…. Новых сеток/каналов/физики нет.

## Выполнение и сохранённые свидетельства

Ровно5entered/completed rawΦ и4entered/completed Jv, native0/Γupdates0.7immutable request records согласованы с caps и completed нижними границами. Metadata всех5Φ образуют sequence1..5 и явно определяют completed как «raw contraction produced», не cone/IO/extraPI success.12rawHDF файлов существуют, каждый151,006,208B<2GiB; это base,4probes,2central,analyticJΓ,4Q. Рецензент проверил только filenames/stat bounds и малые TOML/log/code, HDF содержимое не перечитывал. Whole rawΦ replay ровно0 и все9scalars bitwise — свидетельство actual проверенного adapter, не второй reviewer full-array scan. Девять scalar current/parent значений отдельно побитно перепроверены против parent map064.

Base raw residual тот же исходный: relative4.125583337831235e−6/weighted3.412502855963494e−6, weighted absolute5.555432467297159e−7; это **не** достигнутый trigger1e−9. Replay не превращает этот несошедшийся input в root. Все5probe qcov finite/nonnegative. Ни μ-fit для старого source, ни Φ не заменялись нормировкой/clip; новые current μ входят только в разрешённые5Φ.

## Производная

Для v=Γ, ε=.001/.0005: analytic−FD relative4.700957295750325e−6→1.1752497426512688e−6, weighted4.756628823540782e−6→1.189168985819738e−6. Независимо пересчитанные absolute-error ratios3.99996454/3.99996038 соответствуют ε²; central mutual relative3.5257106431292654e−6. Empirical FD gate действительно проходит. Это сильное проверяемое evidence первого derivative **в одном направлении**, не строгая error enclosure и не проверка всех Arnoldi directions с точностью1e−10.

Все4Jv проходят прежние arithmetic gates. Прямая spectral recomputation relative GR5.4149247428365434e−15/A5.357500578858921e−15; PI≤4.82263e−16; raw/Herm Jv conversion≤4.79648e−16. Validation δNrelative2.7974925797755126e−15; максимум **среди всех4**1.362045357370233e−14 (первый Krylov Jv), а не validation-only значение. χ=.6469060727688077>0; δμ validation−.012081342160135549 в nativeE0; formula−numerator/χ и δNabs/scale независимо совпали. Tiny31assertions PASS, включая independent dense18/45-chain derivative,5 actual native linear channels и V2 serialization/raw-before-guard controls; exit0/12.51s/656196KiB/swap0.

Source initialidentity.method по freeze всё ещё содержит «exact-PD pending». Отдельный уже принятый ALL certificate a703ec4c… из `reviews/global-gamma-exact-pd-result.md` связан Root admission с тем же Γfilecab810…: все337744 сохранённых binary blocks strictPD. Это дополнительное внешнее свидетельство, которое данный reviewer не пересканировал; старую frozenlabel не переписываем. Оно устраняет именно условность PSD-base для scaling probes, но не доказывает full nonlinear convergence/forward accuracy или положительную окрестность заданного радиуса для произвольного Arnoldi direction.

## Малый Krylov oracle

StdlibPython≤128MiB независимо из Hbar.rows(4×3) пересчитал cubic roots3×3, normalized eigenvector cross-products/Ritz residuals, singular values через eigen(BᵀB), least-squares coordinates отдельным partial-pivot solve нормальных уравнений, residual и Bᵀres≈0. NumPy/Julia/newmaps не использовались. Все значения совпали с сохранёнными (relative≤2e−11).

- Ritz=[−.112832345715821,.080069073072304,1.000148055073442]; residuals=[.7819393533,.7581753382,.1168493332]. Разность near1 Ritz от1≈1.48e−4 несопоставима с leakage.11685; full spectral radius/stability/stall mechanism не установлены. Для nonnormal operator нельзя обращаться с этим Ritz как с доказанным eigenvalue.
- B=Ibar−Hbar имеет restricted singular values[1.5727097412,.8786367635,.0869352319]. Сохранённый unit near-null coordinate даёт ||Bz||=.086935231915; это span witness, не полный σmin certificate. H21=.2012,H32=.7488,H43=.7847 — здесь не наблюдается near breakdown/нормализация почти нулевого нового вектора.
- GMRES predicted residual5.467025195624652e−7 против seed5.555432467297161e−7: ratio.9840863385, уменьшение1.591366148%. Это только projected linear prediction; Γ correction не сформирована/не применена. Full saved Q orthogonality и field-level Arnoldi relation reviewer не перечитывал, поэтому full-space заключение опиралось бы на дополнительное field evidence.

Seed Herm conversion relative3.74157e−11/weighted4.32499e−11 сохраняется отдельно от raw residual; weighted absolute2.40272e−17. Это не стирание исходной невязки: actual trigger и reported residual остаются raw.

## Ресурсы и границы решения

Actual log: exit0,82.72s process/81.831379s script, MaxRSS3,444,604KiB=3,527,274,496B,Swap0; external600s/12.125GiB/high11.5/J4/BLAS1/data3/file2GiB соблюдены. Итоговые bytes1,812,332,622; summary data_bytes1,812,332,322 измерены до собственной записи300B, это согласованная разница, не скрытый raw. Этот результат не обещает memory/time следующего более длинного Krylov/fine-grid/fullΦ run.

Машинная проверка `stable_source_jv_result_audit.py`/`stable-source-jv-result-audit.json`:23metadata/math/IO-bound checks PASS, без HDF reads. Приняты derivative consistency и отсутствие сильного улучшения в **этом трёхмерном** span; root/acceleration/механизм стагнации/fine-grid accuracy/эксперимент не приняты. Дальнейшее действие должно выбирать новое различающее свидетельство, не применять сохранённые GMRES coordinates как уже проверенный Γ step.
