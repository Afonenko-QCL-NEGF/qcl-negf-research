# Независимый результат HermΛ warm16

**HermΛ warm16 выполнил16Π и уменьшил raw relative residual в13.02раз до6.83727e−4, но reduced trigger не достигнут.** Initial actual native/KMS/cavity gate прошёл, native exit stage не выполнялся. Наблюдение поддерживает дальнейший ограниченный тест этой trajectory; ни restricted convergence, ни native SCBA/physics acceptance пока не установлены.

## Provenance и исполнение

Primary `equilibrium-gamma-warm-hermitian/attempt01`, imported Coree171/tree742 clean,Julia1.13.3/4threads/BLAS1. Archived maina6199af0/functions3b296973/warmseedd786/checkscc0/warmchecks7262 соответствуют просмотренному freeze; exact source hashes независимо совпали identity. Actualinput SHA256 `6f9813c6445b7b5901539d2835a672ad29cfc05c47a41bdc2d4c30d1a78ca70c`; standardaccepted Manifestda0780fb прежний. Cgroup12.125GiB/high11.5/Swap0/CPU400/FSIZE5GiB записан; runtimeflag былauthorized. Seedmanifest18439792/seed4142e06d/sourcea690 иactualproblemfingerprints прошли повторно.

Root log exit0,165.18s processwall,151.136s harness,peak6853164KiB,swap0. Summaryexecutioncomplete,statuspilot_call_cap,data1976932470B,scientific_acceptancefalse.104small checks прошли доactualmodel. Это штатныйcap, неtimeout/OOM. Builder/failedoriginalwarmprobe costs — отдельныеcampaign costs; этотrun не выдаётся заsinglecost принятогонаучногоresult.

Independent `reviews/hermitian_warm16_audit.py`, JSON/CSV: всеsmallTOML/log/sourcehashes,0.009s,Python17408KiB; никакогоJulia/новогоmap/fullarrayread. HDF checkpoint проверен толькоmetadata:Gamma-first call0/actualidentity, Gamma-last call16/schema/freshAArecipe. Саморавенство fullarrays разныхcheckpointfiles не пересчитывалось.

## Gates, corrections, trajectory

Initialactual-one-map passedtrue,worst5.71293496e−16,retarded_sum2.66495655e−16,sameinputcavityGR0; channelKMS primaries сохранены. Это currentwarmseed/currentμ consumer equivalence, неuniform KMS всехбудущихstates. UnchangedCore native candidate используется какconsumer, ноshared LO/PI primitives не становятсяindependentphysicsoracle.

Вlog18PI_HERMITIAN_CORRECTION entries:16Π localR +2extra initialcontrol local reconstructions. Максимумыabsolute2.21181e−15,unweightedrelative2.62140e−16,weightedrelative2.60910e−16 — намного нижеunchanged1e−10 arithmetic guard. Это измерениеexplicitknownsymmetryrepresentation, неΓclipping. Число18непревышаетcap16Π: они считаютразныеоперации.

Всего16maps/16acceptedattempts: seed1,linear1(map2),AAII14(maps3–16), rejected0/fallback0. МаксимальныйacceptedQRcondition7.3761<1e8,L1=6.62072<10. Guards/finitecurrentμ пройдены, тоесть strictnegativeq obstruction вэтойtrajectory не повторился; minq/allnodemargins отдельно не сохранены и не объявляютсяизмеренными.

| map | raw relative | weighted relative | μ,meV |
|---:|---:|---:|---:|
|1|.008899907754|.008837248111|8.46402562|
|16|.000683726713|.000625088410|8.46009340|

Raw numerator.1542394613→.01183642859,denominator17.33045618→17.31163690. Улучшениеrelative13.0168×,absolute13.0309×,weightedrelative14.1376×; denominatorпочтинеизменен. WorkingHermresidual отличается лишьroundoff; maprepresentationcorrection≈9e−17relative. Следовательноуменьшаетсяrawmaperror, не толькоmixedstep илиnormscale.

Observedrawresidual строгоубывал, но математическаяcontractivity/remainingerror из16точек не доказаны. FinalrΓ вышеtrigger1e−9 примерно6.84e5раз,weighted6.25e5раз; all16triggerfalse, ниодногоpassstreak. Parentbest=last поscalarbest_relative, statusчестноcap16.

Currentμfit — constructivecurrentA constraint: maxrelativefiterror1.974e−13,final2.267e−14,target.045. Егоsmallness/почтипостоянныйμ не являютсяrawKeldyshλ илиequilibriumproof. Finalinput/mapΓ имеютcomputedminima≈−5.14e−19/−5.15e−19,within7.11e−14 backwardallowance; Aminimum≈−2.37e−15 такжеwithinallowance. Noeigenclipping, ноstrictPSD/PD иуниверсальнаяcausality изэтого не следуют. SelectedBigoracle ужепоказал matrixPSDproblem HermΛ-only, отличныйотqtrace; текущийscalartrajectory его неотменяет.

## Сопоставимость и отсутствующая native приёмка

Исходныйcoldattempt02 завершилсяr16=.0203879940 противнынешнего.0006837267, иначальнаяневязкабыла1.9285 против.0088999. Это **не pure warm-start matched control**: одновременно изменилисьinitializer иretardedrealpartrepresentation(HermΛ). Oldwarm безHermΛ остановилсядоΦ, поэтомунепредоставляетtrajectory дляизоляцииwarm-effect. Данныепоказываютуспехкомбинацииwarm+HermΛ; не приписыватьвсе29.82×различиеfinalresidual толькоinitialization. Внутриэтогоrun13×уменьшениеизмерено корректно.

Native summarynot_measuredbecause no_reduced_trigger; materialization/solver/audit events иstarted/completedbounds все0. Значитнетnewnative rK/rΣ/rawλ/candidateλ,per-familyRLG/embeddingidentity/collision/FDR/physicalcurrent audits. Initialnativeconsumerоперации выполнены вprecheck3 и неявляютсяnativeexititerations. FDGreen зановостроилсякаждуюΠ, поэтомуналичиеFDpair не переноситthermalconstraint вnativevalidation заднимчислом. Poisson/U/discretization/experiment не измерены.

## Checkpoint и допустимое продолжение

`Gamma-last.h5` schema campaign133-gamma-checkpoint-v1,pi_calls_total16,recipe **Gamma field only; fresh AA history, not native exact restart**. Fixedwriterсохраняетaccepted inputΓ map16, анеΦ16. `Gamma-best` scalarbest такжеcall16. Приследующемrestart128 должныбытьfrozen parentinput/model/Gamma-file hashes, толькоΓ initializer иemptyAA; initialboundary обязанповторитьmap16 raw/weighted/currentμ втомжеfunctional. ЗагрузкаrawΦ16 вместоΓ16 былабыдругимначальнымполем, незадокументированнымupdate.

Continuation128 — новаяidentity/cap/resource/dataadmission, не ретроактивноеувеличениеэтого16attempt и неexactAAresume. Послеboundary продолжаетсятотжеΦ/localHermΛ operator/currentμ/nativegates; nativevalidationостаетсятолькоусловнойпопрежнемуtrigger. Настоящийнаучныйвыводэтогозавершённогоrun:boundedrawresidualreduction, безconvergenceacceptance.
