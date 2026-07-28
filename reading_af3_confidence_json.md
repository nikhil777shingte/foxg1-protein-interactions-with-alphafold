# Reading an AF3 `summary_confidences_0.json` file

Field-by-field breakdown of what each value in the confidence JSON means, using the
FOXG1-HDAC6 complex as the worked example (a failed/low-confidence prediction, which is
useful for seeing what the "bad" numbers look like).

Companion to `script_measurement_summary.md` and `reading_lysine_contact_output.md`.

## Example file

```json
{
 "chain_iptm": [0.19, 0.19],
 "chain_pair_iptm": [[0.25, 0.19], [0.19, 0.64]],
 "chain_pair_pae_min": [[0.76, 19.88], [20.14, 0.76]],
 "chain_ptm": [0.25, 0.64],
 "fraction_disordered": 0.44,
 "has_clash": 0.0,
 "iptm": 0.19,
 "num_recycles": 10.0,
 "ptm": 0.52,
 "ranking_score": 0.47
}
```

## Field-by-field

| Field | Value | What it means |
|---|---|---|
| `chain_iptm` | `[0.19, 0.19]` | Each chain's overall interface confidence. Equal for both here — normal for a 2-chain complex, since it's the same interface viewed from either side. |
| `chain_pair_iptm` | `[[0.25, 0.19], [0.19, 0.64]]` | A small 2×2 grid. Diagonal = each chain scored against *itself* (0.25 = FOXG1 alone, 0.64 = HDAC6 alone — matches `chain_ptm`). Off-diagonal (0.19, both directions) = the actual FOXG1↔HDAC6 interface score. |
| `chain_pair_pae_min` | `[[0.76, 19.88], [20.14, 0.76]]` | Same grid, but positional error in Å. Diagonal ≈ 0.76 Å (trivially low — a chain compared to itself). Off-diagonal ≈ 20 Å — the best-case interface positioning error. |
| `chain_ptm` | `[0.25, 0.64]` | How well each chain folds *on its own*, ignoring the other chain entirely. FOXG1 = 0.25 (poor — matches its disorder), HDAC6 = 0.64 (decent). |
| `fraction_disordered` | `0.44` | 44% of all residues in the complex (both chains combined) are predicted to have no fixed structure — the floppy loops visible in the structure viewer. |
| `has_clash` | `0.0` | A coarse whole-chain clash check. 0 = no gross overlap detected. Doesn't catch small individual atom-pair clashes (e.g. sub-2 Å contacts) — it's a blunt instrument, not a precise one. |
| `iptm` | `0.19` | The single headline interface score for the whole complex. Same number as the off-diagonal in `chain_pair_iptm`, since there's only one interface to score in a 2-chain complex. |
| `num_recycles` | `10.0` | Not a confidence score — a setting. How many internal refinement passes AF3 ran before finalizing this structure. |
| `ptm` | `0.52` | Overall structural accuracy for the *whole complex as one blob* (both chains and the interface, mixed together). Less specific than ipTM for judging the interface alone. |
| `ranking_score` | `0.47` | AF3's composite score combining ipTM, pTM, clashes, and disorder into one number — this is what decided which of the 5 generated models got called `model_0` (the top-ranked pick). |

## Which two numbers actually matter

For judging whether a predicted interface is trustworthy, everything above reduces to two
numbers:

- **`iptm`** — the headline interface confidence (0.19 here — well below the 0.6 threshold
  AF3 uses to separate "gray zone" from "likely failed")
- **The off-diagonal value in `chain_pair_pae_min`** — the interface positioning error in Å
  (~20 Å here — AF3 doesn't know where FOXG1 sits relative to HDAC6)

Everything else in the file either restates those two in a different shape (`chain_iptm`,
`chain_pair_iptm`), explains *why* they're low (`chain_ptm`, `fraction_disordered`), or is
an unrelated technical setting (`num_recycles`, `has_clash`, `ranking_score`).

## For comparison: a better-scoring complex

```json
{
 "chain_pair_iptm": [[0.26, 0.40], [0.40, 0.76]],
 "chain_pair_pae_min": [[0.76, 8.44], [9.65, 0.76]],
 "iptm": 0.40
}
```

(FOXG1-HDAC1, abbreviated.) Off-diagonal ipTM of 0.40 and PAE-min of ~8-9 Å — still under
the 0.6 confidence threshold, but a meaningfully more coherent prediction than the HDAC6
case above. See `results_caveats_section.md` for how this distinction was used when
writing up results for researchers.
