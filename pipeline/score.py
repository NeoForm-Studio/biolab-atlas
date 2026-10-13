"""Extremophile scores (0-100) and how well they agree with lab measurements (BacDive).

heat: position of the estimated growth temperature among all profiled species.  cold: 100 - heat, so it separates
cold- from heat-adapted genomes, not psychrophiles from ordinary mesophiles.  salt: position of proteome acidity,
raised to 90 when the full ectoine pathway is present.  radiation: share of the 6 Deinococcus markers found.
Heat and cold are checked against BacDive temperatures with ROC AUC (0.5 = chance, 1 = perfect). Most labels are
culture temperatures, and the species were chosen by hand, so the check mostly tests extremes.
"""

THERMOPHILE_C, PSYCHROPHILE_C = 50, 20  # reference-temperature cutoffs for the validation labels


def percentile(values, v):
    """Mid-rank percentile of v among values (0-100)."""
    below = sum(x < v for x in values)
    equal = sum(x == v for x in values)
    return round(100 * (below + equal / 2) / len(values))


def auc(positives, negatives):
    """Mann-Whitney ROC AUC: chance that a random positive outscores a random negative (ties count half)."""
    if not positives or not negatives:
        return None
    wins = sum((p > n) + (p == n) / 2 for p in positives for n in negatives)
    return round(wins / (len(positives) * len(negatives)), 2)


def reference_temp(row):
    """BacDive optimum midpoint if known, otherwise the midpoint of the temperatures it was cultured at."""
    for key in ("measured_optimum_c", "measured_growth_c"):
        r = row.get(key)
        if r:
            return (r[0] + r[1]) / 2
    return None


def upper_temp(row):
    """Upper end of the temperatures BacDive reports the strain growing at, else of its optimum."""
    for key in ("measured_growth_c", "measured_optimum_c"):
        r = row.get(key)
        if r:
            return r[1]
    return None


def found_fraction(text):
    """'5/6' -> 5/6; anything else -> None."""
    try:
        a, b = text.split("/")
        return int(a) / int(b)
    except (AttributeError, ValueError, ZeroDivisionError):
        return None


def score_library(rows):
    """Add a 'scores' dict to every row that has a proteome; return the validation summary."""
    ogt = [r["predicted_ogt_c"] for r in rows if r.get("predicted_ogt_c") is not None]
    acid = [r["acidic_excess"] for r in rows if r.get("acidic_excess") is not None]
    for r in rows:
        if r.get("predicted_ogt_c") is None:
            r["scores"] = None
            continue
        heat = percentile(ogt, r["predicted_ogt_c"])
        salt_out = (found_fraction(r["markers"].get("salt")) or 0) == 1  # full ectoine pathway
        r["scores"] = {
            "heat": heat,
            "cold": 100 - heat,
            "salt": max(percentile(acid, r["acidic_excess"]), 90 if salt_out else 0),
            "radiation": round(100 * (found_fraction(r["markers"].get("radiation")) or 0)),
        }

    labelled = [(r, reference_temp(r)) for r in rows if r.get("scores") and reference_temp(r) is not None]
    hot = [r["scores"]["heat"] for r, t in labelled if t >= THERMOPHILE_C]
    not_hot = [r["scores"]["heat"] for r, t in labelled if t < THERMOPHILE_C]
    cold = [r["scores"]["cold"] for r, _ in labelled if upper_temp(r) <= PSYCHROPHILE_C]
    not_cold = [r["scores"]["cold"] for r, _ in labelled if upper_temp(r) > PSYCHROPHILE_C]
    return {
        "heat": {"auc": auc(hot, not_hot), "positives": len(hot), "negatives": len(not_hot),
                 "label": f"BacDive reference temperature >= {THERMOPHILE_C} C"},
        "cold": {"auc": auc(cold, not_cold), "positives": len(cold), "negatives": len(not_cold),
                 "label": f"upper end of BacDive temperature range <= {PSYCHROPHILE_C} C"},
        "salt": {"auc": None, "note": "not validated: BacDive salt records are too sparse and mix tolerance with need"},
        "radiation": {"auc": None, "note": "not validated: no lab radiation measurements in BacDive; "
                                           "markers only detect the Deinococcus strategy"},
    }
