# Future raw-response control V2: fresh-base anchor и prepaired access

Только уточнение future proposal по независимой review
`reviews/noEE-frozen-response-and-linewidth-review.md`; V1
`next-noEE-frozen-response.md` SHA68bf83f816bae77ff02a8b9ce4d5793a3067b49a7da3fba93298e8185252c56b
не изменён. Это не approval implementation/actual. Julia/Pkg/H5/опытов нет.

Исправление обязательное: frozen arm использует **R* от fresh unperturbed
full map**, а не stored input R. Сначала из source G0 вычислить свежие full
Σ*R/L/G, side-cavity retarded coefficients/embeddings и output centralGR*.
Сохранить raw base n*,p* **до** paired normalization. Обе frozen applications
(G0 и perturbedG0) используют именно эти local/side/central retarded coefficients.

До интерпретации δ-response требуется gate: frozen-base **raw** n,p whole-field
воспроизводят fresh-full-base raw n*,p* при одном объявленном machine arithmetic
допуске (absolute/scale/weighted и bitwise flags отдельно). Если replay не
проходит, control unavailable; обычный base-cycle drift нельзя назвать
occupation→R feedback. Лишь после replay сравнивать full-perturbed response
с frozen-perturbed response, вычитая общий fresh-full-base raw output.

`_current_green_production` на e171 уже возвращает **normalized** Green.
Нельзя использовать его GL/GG как raw. До actual нужны конкретный reviewed
raw-capture adapter или точное reconstruction из **тех же** outputGR* и total
Σ*L/G с неизменным product order; full/frozen base replay должен проверить
также этот raw-access contract. Capture до `_paired_number_coefficients`
предпочтителен; fresh source arrays не заменять позже stored/mixed families.
Сохранить raw Ns/Ps, a,c, `Nnew_s=aNraw_s+cPraw_s` и paired correction отдельно.

One perturbation/globalN/sourceB preservation и cap≤2full maps+2frozen-LG
applications остаются. Fine48001/25µeV потребует **новой source identity и
graph**: ideal1468/40 даётgcd4; parity тогда coarse contrast, не доказанные
ровно2components. Budget/доступ к raw пока не реализованы и не разрешены этой
запиской. Такая постановка локализует mixed/full response, не доказывает
причину stall, spectral radius или physical sector conservation.
