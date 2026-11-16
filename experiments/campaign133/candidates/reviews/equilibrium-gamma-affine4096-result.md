# Независимый результат affine4096 pilot32

**Более высокий affine cap позволил принимать новые AA steps, но не дал сходимости.** За32TOTALΠ raw rΓ уменьшился7.19231e−6→5.81009e−6, best5.80635e−6 наcall30. Семь QR-admissible AA proposals остановил PSD guard. Native exit не достигнут. Это полезный отрицательный различающий результат: L1cap10 ограничивал часть предложений, но простое увеличениеcap не устранило все препятствия.

## Provenance, replay, стоимость

Primary `equilibrium-gamma-affine4096/attempt01`: frozen functions cf923bfb/restartb020d599/main876a8c3e соответствуют independentlychecked sourceidentity; actualinput и всеarchivedcodehashes совпалиidentity. SourceΓlast128b30faef1/finalparent57f2/inputfddd/model47682/pinnedmap128 прошли loader. Coree171/tree742,Julia1.13.3/4threads/BLAS1/caps12.125GiB/high11.5/swap0/CPU400 unchanged. Parent не подменён live/current selector.

Boundary replay passed: whole rawΦ difference0 (weighted/unweighted/maxblock), все9scalar fields exactlyparentcall128. Newcall1raw7.192305218634285e−6; currentμ занововычислен иexactparent. Initialnativeconsumer worst5.72106e−16; это currentinputcontrol, неuniform futureKMS/nativepass. EmptyAA/streakreset намерен, поэтому неexactAAresume.32calls содержатоднуreplayedboundary и31дальнейшуюevaluation.

Runexit0,270.27s processwall,260.832s harness,peak6582376KiB,swap0,data3034020864B. Summarypilot_call_cap/scientific_acceptancefalse, nativeevents/start/completionbounds0,reasonno_reduced_trigger. Новые24tests passed; oldphysical104/restart15 evidence reusedявно. Это неtimeout,OOM,IODiskrefusal илиacceptedscientificresult.

Independent smallaudit: `reviews/gamma_guard_results_audit.py`, `gamma-guard-results-audit.json`,CSV отдельнодляparent128/new32. ТолькоTOML/log/sourcebytes, noJulia/maps/largearrays; `/usr/bin/time` repeat.07s,peak19476KiB. Initialru_maxrssаномалиятогожеtoolcontext880128KiB не использована каквалиднаяcostmeasurement; численныепоказатели повторились. Финальныйparent128отдельноразобранв`equilibrium-gamma-cont128-result.md`.

## Какие proposals действительно приняты

32mapfiles/39attemptrecords: acceptedseed1,linear8,AAII16,convexfallback7; rejectedAA7. Всеrejectionsпроизошли **до Π** наΓcone, поэтому32actualΠbudgetнепревышен, failed evaluatedΠнет. Nativeprecheck3 — отдельныйсчётчикопераций, неnativeexititerations.

Из16acceptedAA пять имелиL1>10:32.04(call8),119.42(13),20.74(19),19.55(25),396.06(26). Этиproposals старыйcap10непропустилбы; фактический механизмизмененияalgorithmadmission подтверждён. ОстальныеacceptedAAимелиL1≤10; finalL1≈1.365,rank3,QRcondition7.246.

СемьотклонённыхAAпрошлиQR/L1≤4096, ноимелиnegativeΓeigenvalue:

| передследующимmap, budgetcall | L1 | λmin proposal | node(e,k) |
|---:|---:|---:|---|
|5|2389.17|−3.056e−13|(5090,1)|
|8|3530.37|−2.320e−13|(5090,1)|
|13|1592.86|−1.704e−13|(5090,1)|
|16|2906.38|−3.342e−11|(9771,1)|
|19|1595.01|−2.568e−11|(9771,1)|
|22|585.78|−2.245e−12|(9771,1)|
|26|355.56|−5.485e−12|(9771,1)|

Declaredallowance7.10543e−14: нарушения≈2.4–470разпревышают его. Этонеоснованиеclip'итьΓ илиослабитьcriterion. FloatLAPACK/domainguardотклонилвычисленныеproposals; истиннаяпричинаихнаправления — nonlinear geometry, amplifiedmap/inputnoise илисочетание — этимscalarlogнеразделена. Allrejectionsсохранилиrawtrial; одинconvexfallbackкаждыйразaccepted,historyreset.

Suffixtelemetry надо считатьбезdoublecount: acceptedfallbackповторяетinfo исходногоrejectedproposal. Послеудалениятакихдубликатов23QR-admittedproposals,6rejectedlongersuffixesпоL1(≈8752–45643),ниодногоrank/conditionrejectвэтойtrajectory.23=16actualacceptedAA+7PSDreject. Вparent128105`linear_rank_reset` не доказываютrankdeficiency: староеимяобозначалоfallbackпослевсехguards, detailedsuffixcausesнет.

## Raw progress и пределы интерпретации

| call | raw relative | weighted relative |
|---:|---:|---:|
|1|7.19230522e−6|5.95754891e−6|
|16|7.07813857e−6|5.87687753e−6|
|30(best)|5.80634820e−6|4.83409104e−6|
|32(last)|5.81008885e−6|4.83725590e−6|

Finalimprovement1.23790× поraw и1.23160×weighted, тоестьrawуменьшился≈19.2%. Bestчутьлучшеlast; есть10upwardrawtransitions, nonmonotoneAAэторазрешает. Ниplateau/короткаяoscillation, ниуменьшениерезидуала не доказываютcontractionconstant илиdistance-to-root. Finalrawещё5810×trigger1e−9, all32triggerfalse.

Raw numerator.000124509567→.000100581359,scale17.3114965→17.3115010: improvementнеотdenominator. Currentμ8.459613669→8.459612207meV; maxdensityfiterror1.67305e−13 — constructiveconstraint, неnativeλ.34HermΛcorrectionlogs=32map+2initialcontrols,maxrelative2.59789e−16/weighted2.59600e−16. RawΦrepresentationrelative≤9.43373e−17; самэтотknownHermitianconversionнеобъясняетglobalresidual≈6e−6. НоfullΠaccuracy несертифицирована FFTwitness/μfit поотдельности, особеннопослеlargecoefficients.

Прямойmatchedcounterfactual32calls сcap10 **изтогожеΓ128 иfreshhistory не запускался**. Поэтому19.2%progressнельзяцеликомприсвоитьcap4096 противextraevaluations/historyreset. Доказанболееузкийmechanism: acceptednewL1>10 иреальныйconeobstructionпоследальнейшегоослабленияL1guard. Blindдальнейшееувеличениеcapэтимнеобосновано.

## Следующий дешёвый различающий контроль

Предпочтительнееещёнегонятьдлинныйcontinuation, асделатьотдельныйfixed-input diagnostic:

1. Изсохранённыхrawrejectedtrial выбратьworste9771/k1 иnear-budgete5090/k1: независимыймалыйHermitianeigen/Rayleighaudit соscale/residual/increasedprecision отличитLAPACKsignuncertainty отrobustout-of-coneвычисленногоproposal. ОнсамнеустановитточностьвсейΠ.
2. ЗатемнаfrozenacceptedΓ (новаяidentity) получитьодинAAproposal снебольшимfixedbudget иизмеритьмаксимальныйadmissiblestepθ вдольΓ+θ(ΓAA−Γ), **безeigen-clipping**. ТолькооднаfreshΠнавыбраннойfeasibleточке, плюссопоставимыйlinearα.2point, дастразличение: полезноеAAнаправлениеблокируетсяconeboundary илидажеfeasibleнаправление неменяетrawdefect. Приоченьмаломθконкретнопоказанgeometrybottleneck; приO(1)θснеулучшающимсяrawresidualнужноисследоватьslow/nonlinear modes, приусловииотдельнойоценкиΠnoise. Самоnegativeproposal/plateau не доказываетnoncontractive Jacobian.

Этоproposal диагностическогоалгоритма, неразрешениеегоавтоматическогозапуска и невыборновыхscientifictolerances. Πnoise желательноизмерятьнаодномΓсindependentarithmetic/precisioncontrol, ане издетерминированногоexactboundaryreplay; последнийдоказываетповторяемость, неaccuracy.

## Неизмеренное и checkpoint lineage

Нетnative rK/rΣ/rawλ/candidateλ/individualRLG/embedding/FDR/collision/current audits; physicalPSDconesимеютabsoluteallowance, строгоеPDнеустановлено. Poisson/grid/window/basis/k/experiment outcomesнеизмерены. СамinitialKMSconsumerpass ихнеподменяет.

Gamma-last — inputmap32. **Gamma-best — inputmap30**, хотя genericmetadata`pi_calls_total=32` фиксируетbudgetнаконецзаписи, анеiterationbestfield. Будущийrestartbestнельзяreplay'итьпротивmap32/terminalrawΦ32; нужноявноеbest-call/sourceevidence. CurrentΓlast-onlycontractэтойошибкинеимеет. Scientificoutcome:boundedalgorithmprobeснебольшимулучшениеминеснятымcone/noise/modeвопросом, безconvergenceacceptance.
