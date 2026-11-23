| Case | D | tables (s, MB) | ours: s per matvec of 4 | Gmadd/s | KSHELL: s per matvec of 4 | CoSMo stored, 1 vector |
| --- | --- | --- | --- | --- | --- | --- |
| 24Mg_sd | 28,503 | 0.351 s, 3.5 MB | 0.0036 | 8.29 | 0.00949 | - |
| 28Si_sd | 93,710 | 0.061 s, 6.7 MB | 0.0136 | 9.84 | 0.02553 | - |
| 48Cr_kb3g | 1,963,461 | 0.75 s, 101.9 MB | 0.5804 | 9.68 | 0.70572 | - |
| 24Mg_1hw- | 978,948 | 1.382 s, 130.5 MB | 0.3015 | 9.3 | 0.55611 | 0.073 |
| 28Si_1hw- | 3,834,452 | 3.2 s, 291.0 MB | 1.4017 | 10.13 | 2.18438 | 0.365 |

| Case | ours: matvecs, solve s | KSHELL: matvecs, solve s | max abs dE (MeV) |
| --- | --- | --- | --- |
| 28Si_sd | 500, 2.44 | 352, 3.06 | 1.8e-08 |
| 48Cr_kb3g | 676, 119.79 | 412, 84.53 | 4.0e-08 |
