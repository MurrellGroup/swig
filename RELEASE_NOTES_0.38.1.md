# Swig 0.38.1

## Cross-regime AER-R R-optimized correction

- Retunes the opt-in AER-R `r_optimized` profile against both supplied 20,000-record KIMDB macaque simulations rather than the lower-SHM regime alone.
- Restores the V/D/J mismatch cost from −4 to −3 so high-SHM templated edge tracts remain favorable.
- Uses V scoring `+2/−3/−9/−1`; the −9 gap-open cost retains the three original V-stop regression fixes without the observed 42-nt complex-indel V-end failure at −11/−13.
- Uses D scoring `+2/−3/−13/−1`, retains the 5-nt exact-run floor, and changes the D-state cost from 10→8 to 12→10. The two-point relaxation now requires raw D score ≥20 or exact support from at least two distinct locus-matched D template sequences.
- Uses J scoring `+2/−3/−17/−2`. The same-span one-point D ambiguity window and the AER-R candidate-search/joint-partition algorithm are unchanged.
- Adds the complex-indel IgG V-end case to the deterministic WASM regression fixture.

## Accuracy summary

Against Swig 0.37.6 plain AER, the revised profile improves full-data V/D/J ambiguity-aware Brier loss by 5.60%, 8.69%, and 6.32% on the lower-SHM simulation and by 3.76%, 10.73%, and 5.02% on the IgG-like simulation. V-end, D-start, D-end, and J-start mean absolute errors improve on both complete simulations. These are simulator-specific development results, not external biological validation.
