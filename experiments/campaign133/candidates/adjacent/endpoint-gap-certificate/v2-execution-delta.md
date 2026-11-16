# Endpoint V2: execution-only BLAS verification

V1 `endpoint_certificate.jl` SHAe031d0c7…40d87e3 и весь frame сохранены в
`unexecuted-v1/` вместе с исходным SHA256SUMS и archive-local manifest.

Setter `BLAS.set_num_threads(1)` уже присутствовал в V1 **до** source_frame/
interval work. V2 добавляет только getter до setter, getter после, assertion
`after==1` и запись before/after/Julia thread count в summary. Environment
`OPENBLAS_NUM_THREADS=1` не считается свидетельством фактического thread count.
На source gate failure до output metadata состояние остаётся в root log;
никакого BLAS/physics acceptance по одному environment value нет.

Все interval formulas/constants/kernel tail/Gershgorin/geometry/areas,
frozen sources/PD/checkpoint/hH и caps2GiB/J1/BLAS1/CPU1/120s неизменны.
Автор V2 Julia/Pkg/H5 не запускал. Independent reviewer может продолжать
numerical review V1 и отдельно проверить эту execution-only дельту.
