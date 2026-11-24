# Data dictionary

Metadata for the three top-level CSVs in this folder. All files use comma
separators, a header row, and `NA` as the missing-value code. Column names below
are shown as written on disk; note that `read.csv()` replaces spaces with dots
(e.g. `Number of visits` -> `Number.of.visits`, `Other wild bee` ->
`Other.wild.bee`) and appends a dot to the trailing space in `Crop `.

Licence: [CC BY 4.0](../LICENSE-data).

---

## `VisitationRates_Data.csv`

Observational pollinator-visitation data compiled from five published studies,
used to fit the saturating seed-set model per crop.

4306 rows x 6 columns. No missing values.

| Column | Description | Values / units |
| --- | --- | --- |
| `Study` | Source study identifier | `Chirila_etal_2026` (lima bean), `Kandori_etal_2024` (strawberry), `Kendall_etal_2021` (blueberries), `Garratt_etal_2016` (apple), `DiTrani_etal_2024` (watermelon) |
| `Crop` | Crop common name (dots for spaces) | `Lima.bean`, `Strawberry`, `Southern.Highbush.Blueberry`, `Rabbiteye.Blueberry`, `Northern.Highbush.Blueberry`, `Apple`, `Watermelon` |
| `Pollinator` | Pollinator taxon or treatment | 24 values. Taxa use dots for spaces (`Apis.mellifera`, `Tetragonula.carbonaria`, ...). `Control` / `Unpollinated.Control` = no visits; `mix` / `pollinator.mix` / `Unrestricted visits` = open-pollination treatments |
| `Number of visits` | Number of pollinator visits applied | integer, 0-30 |
| `PollinationMetric` | Type of response measured | `Seed.set`, `Fruit.set`, `Fruit.weight`, `SeedNumber` |
| `PollinationValue` | Response value | units depend on `PollinationMetric` (see below) |

Units of `PollinationValue` by metric:

- `Seed.set`, `SeedNumber`: number of seeds (range 0-470).
- `Fruit.set`: binary, 0 or 1.
- `Fruit.weight`: fruit weight in grams (range 0-3.46). Converted to seed count
  in `scripts/Final_analysis.R` via `WtoS = (w - 0.09) / 0.06`.

## `VisitationRates_ParameterData.csv`

Literature-derived crop parameters: selfing rate (`a`), maximum seed set /
ovule number (`b`), and single-visit efficiency per pollinator group, for 27
crop rows (26 crops; the `chickpeas, lentils` row is dropped in the analysis).

27 rows x 17 columns.

| Column | Description | Values / units |
| --- | --- | --- |
| `Crop ` | Crop common name (trailing space) | text |
| `scientific name` | Latin binomial | text |
| `family` | Plant family | text |
| `a_siopa` | Selfing rate from Sipoa | proportion (0-1); `NA` = unknown |
| `Cultivar_a` | Selfing rate for the cultivar (`a` estimate) | text, `NA` = unknown |
| `a_seed_set` | Seed set under selfing (no visitation) | seeds |
| `Cultivar_b` | Maximum seed set for the cultivar (`b` estimate) | text, `NA` = unknown |
| `b_ovules` | Maximum seed set / ovule number | seeds |
| `Cultivar` | Single-visit efficiency for this cultivar | text, `NA` = unknown |
| `Honeybee` | Single-visit efficiency, for honeybees | seeds per visit, `NA` = no data |
| `Bumblebee` | Single-visit efficiency, for bumblebees | seeds per visit, `NA` = no data |
| `Other wild bee` | Single-visit efficiency, for other wild bees | seeds per visit, `NA` = no data |
| `Hoverfly` | Single-visit efficiency, for hoverflies | seeds per visit, `NA` = no data |
| `Other` | Single-visit efficiency, for other visitors | seeds per visit, `NA` = no data |
| `Reference DOI_a` | Source for `a` | DOI or URL, `NA` = none |
| `Reference DOI_b` | Source for `b` | DOI or URL, `NA` = none |
| `Reference DOI_c` | Source for single-visit efficiency | DOI or URL; several entries are "guesstimate by Nacho" = expert estimate, not literature |

## `parameters_real_crops.csv`

Model parameters fitted to the observational data (`VisitationRates_Data.csv`),
comparing the fitted values (`a`, `b`, `c`, `OVR`) with values predicted from
literature parameters (`a_model`, `b_model`, `c_model`, `OVR_model`).

7 rows x 13 columns. `a_model`, `b_model`, `c_model`, `OVR_model` are `NA` for
the two blueberry varieties (`NHighbushBlueberry`, `Rabbiteye.Blueberry`).

| Column | Description | Values / units |
| --- | --- | --- |
| `` (empty) | Row index written with a leading comma | integer |
| `X` | Short row label | text |
| `crop` | Crop common name | `apple`, `strawberry`, `lima bean`, `SHighbushBlueberry`, `NHighbushBlueberry`, `Rabbiteye.Blueberry`, `Watermelon` |
| `D_R2` | Difference in R2 (linear/quadratic minus mechanistic model) | R2 units |
| `D_AIC` | Difference in AIC (linear/quadratic minus mechanistic model) | AIC units |
| `a` | Fitted selfing parameter | % of `b` |
| `b` | Fitted maximum seed/fruit set | seeds |
| `c` | Fitted diminishing-returns parameter | unitless |
| `OVR` | Optimal visitation rate (fitted) | visits |
| `a_model` | `a` predicted from literature parameters | % of `b` |
| `b_model` | `b` predicted from literature parameters | seeds |
| `c_model` | `c` predicted from literature parameters | unitless |
| `OVR_model` | Optimal visitation rate (literature) | visits |

---

## Other data

`Watermelon_data/` holds the raw watermelon records and the script
(`Extract_watermelon.R`) that created the watermelon rows in
`VisitationRates_Data.csv`.
