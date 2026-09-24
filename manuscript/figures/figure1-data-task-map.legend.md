# Figure 1 legend

| Corpus | Task | Records | Use | License | Note |
|---|---|---:|---|---|---|
| Temporal absolute-kinetics pool | Absolute kinetic constant | 192 | never scored | CC0 / CC-BY-4.0 per source | Below its 300-record and 30-family gates; predictions are not permitted. |
| EnzEngDB v1 | Relative engineering fitness | 6,423 | zero-shot external | CC-BY-4.0 | 51 homology-cold campaigns; macro-averaged Spearman 0.071. |
| Lunzer IMDH landscape | Absolute kinetic constant | 512 | zero-shot external | CC0-1.0 | 67-69% identity to training homologs: a mutation-sensitivity test, not homology-cold. |
| SABIO-RK Km rows | Absolute kinetic constant | 0 | blocked | not established | Permission request drafted but never sent. |
| Frozen homology-cold test | Absolute kinetic constant | 1,646 | internal test, already observed | derived from training corpus | All headline RMSE numbers come from here. Not independent: it has been scored. |
| CataPro retrained | Absolute kinetic constant | 1,646 | same-split comparator | code MIT; data unresolved | 1,127 of 1,646 test pairs overlap CataPro's published ten-fold data. |
| UniKP retrained (Mode B) | Absolute kinetic constant | 1,646 | same-split comparator | code MIT; data unresolved | Retrained on the frozen training partition, so it shares the test. |
| DLKcat retrained (Mode B) | Absolute kinetic constant | 1,646 | same-split comparator | GPL-3.0-only | Retrained on the frozen training partition. |
| UniKP/DLKcat kcat corpus | Absolute kinetic constant | 16,838 | training | unknown | The training corpus itself; by definition not independent. |

Read the panel as: no corpus sits in the top-left where an independent, homology-cold, never-scored benchmark would belong. The internal test has been scored, the comparators share it, EnzEngDB carries a different endpoint, and the temporal pool has never been scored but is below its own gate.