NOTE (2026-10-06): how to read these D7 results.

RAW-DISTANCE power (mode = raw): interpretable. Diluting the real Normal->Adenoma difference by swapping
samples between groups lowers the test z-score steadily (100% -> 25%), as intended.

SHAPE-ONLY power (mode = shape-only): NOT interpretable as power to detect a smaller effect. Swapping
samples makes each group a mixture of normal and adenoma samples, i.e. two sub-clouds. The shape-only test
detects that two-cluster structure, so the median z often RISES after swapping (e.g. GSE41258 adjusted:
1.51 at 100% -> 4.04 at 75%; GSE37364 unadjusted: 1.63 -> 3.26). The diluted groups contain an artificial
shape difference that the real data do not have. The shape-only power figures and the shape-only
classifications in classification.csv must not be used.

Power at 100% (no swaps) is simply the detection rate of the real Normal vs Adenoma comparison, re-run 50
times with 99 permutations; it shows how reliably that real difference is detected.
