# Swig 0.26.1

## UCA default changes

- Changed the phylogenetic-UCA template-leakage default from `0.003` to `0`. For numerical stability, a selected value of zero continues to be evaluated at the existing `1e-9` arithmetic floor.
- Changed the default phylogenetic-UCA inference route from conditional maximum likelihood to continuous Gibbs/MH.
- Conditional ML and explicit grid marginalization remain available without methodological changes.
- Existing saved sessions retain their saved UCA settings. New analyses and **Reset all UCA settings to defaults** use zero leakage and Gibbs/MH.
