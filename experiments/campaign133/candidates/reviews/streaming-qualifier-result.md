# Streaming spectral qualifier: независимая проверка результата

**Решение: arithmetic equivalence подтверждена для всей зафиксированной Γlast64 на coarse grid.** Streaming helper может использоваться как квалифицированный кандидат этой spectral-фазы в следующем отдельно ограниченном whole-map опыте с сохранением source/native/domain guards. Это не принятие SCBA, не доказательство fine-grid adequacy и не универсальная bitwise-теорема для будущих Γ.

Артефакт `streaming-spectral/qualified-attempt01/qualifier-summary.toml` SHAc48d6e4ff1d1259356933d44702d74db353615eb7c78e57adc1e5cc4b8c65562. Frozen wrapper72cd695c…/helper8df7885c… неизменны; actual runtime input47d0b7f655ea7af97175f271f8bd772911abff2e81aa9fb652c83bdb49c85abc. Protocol после документальной поправки0a15d3b0… скопирован в actual output и совпадает по hash. Правильная lineage: Γlast64 — **вход Γ map064**, тот же TerminalGamma; это не новая Γ после последнего map. Старый bedec protocol сохранён какv1.

Источник: completed stableA-cont64, Γfilecab81029…/model700ae176…, Coree171/tree742, Julia1.13.3/J4/BLAS1; immutable source fields повторно hash-проверены до/после, Γarray digestf5a3778b… не изменился. Рецензент сравнил actual metadata/code/input hashes с Root digest JSON; большие sourceΓ/model повторно не читал. No Φ, μ-fit, FD/mixing/Poisson: структурно и по журналу ровно2 responses started/completed, Phi0, mu_fits0, mixing0.

## Свидетельство эквивалентности

Все5 **полных** массивов GR/A/q/condition/scale имеют difference absolute/relative/weighted/maxblock=0 и bitwise_equal=true. Все10 selected split quantities также имеют0/bitwise=true: Γ, Λraw, ΛHpreassembly, localR, totalR, GR, Acov, Asub, Kphysical и Γeffective. Original1e−10 и planned1e−12 пройдены во всех15 сравнениях. Fullarray утверждение опирается на ранее статически проверенный all-entry comparator и этот неизменный actual report; reviewer полные response arrays не перегружал.

Независимая малая проверка посредством binary H5dump перечитала **пять сохранённых complex-полей** Γ/Λraw/ΛHpreassembly/GR/Acov из обоих selected H5. Все их Float64 bytes совпали; отдельно восстановленные unweighted/wE·wk нормы совпали с report≤2e−13. Это44 blocks: E indices[382,525,8040,30704]×11k, Nb5. Пример: ||ΛH||=.05520749344368878; ||GR||=86.50690986797645; ||Acov||=1.2058074876412704e−11. Самые слабые источники сохранены раздельно: ||selectedΓ||=4.771386406489225e−15, а не извлечены вычитанием из localR. Объём selected-response scope не следует смешивать с полной Γ/model или новыми44 Big доказательствами.

Обе PI roundoff=3.4982122613568972e−16 (без нового критерия). h correction relative9.91225e−17 и HermΛ correction2.50277e−16/weighted2.50584e−16 одинаковы. qcov обеих ветвей положителен: min1.3284020439152628e−43, capacity28.643308967036983; tiny qsub negative tail64 остаётся прежним диагностическим контрастом, не удалённым входом. Нового fitμ или плотности не измерялось.

Scalar diagnostic reductions по блокам дают последние-bit различия относительно old full norms, что было объявлено: например absolute Acov−Asub old8.432277177102614e−12 vsstream8.432277177118324e−12. Это различие порядка редукции при **побитно одинаковых** A/GR, а не изменение spectrum. Same physical-source norm differences аналогичны.

## Стоимость и границы

`qualifier-attempt01.log`: exit0; process25.10s, script20.174175s; MaxRSS3,082,304KiB=3,156,279,296B≈2.9395GiB, Swaps0. Summary Sys.maxrss совпадает с process log точно. Actual admission hard12.125/high11.5GiB/300s/J4/BLAS1/CPU≤400% и file/data32MiB сохранён; actual output1,457,053B. Старые7.5115s/новые3.9595s — одна first-call пара с разным JIT/order/full diagnostics/H5 I/O; **speed ratio/честный benchmark не выводится**.

4F=540,390,400B относится к live fields собственного stream sweep. Здесь сохранены oldGR/A, поэтому paired sweep6F=810,585,600B плюс планы/cache/model/PI/κ/scale/q/JIT/allocator; это не RSS bound. Peak RSS относится ко всему mixed old+new process, не изолированному stream. Native full guard не обходился: source cache estimate12,936,559,056B≤numerical13,019,119,616B. Fine2x/4x, AA history и LO/Φ memory/time этим опытом не измерены.

Root tests независимо просмотрены по журналам: stream62/62 PASS (18.64s process/711888KiB) и comparator6/6 PASS (3.13s/461560KiB), оба exit0/swap0. Рецензент не запускал Julia/Pkg/научный расчёт. Малый stdlibPython+selected H5dump audit≤128MiB: `streaming_qualifier_result_audit.py` и `streaming-qualifier-result-audit.json`; все13 metadata/comparison checks и process log checks пройдены.

Следующий различающий опыт — whole-source map с этой заменой spectral-phase и прежними rawΦ/residual/native gates, отдельной identity/admission и прямым сравнением; перенос coarse битовой эквивалентности на другой grid/Γ/history заранее не принимается.
