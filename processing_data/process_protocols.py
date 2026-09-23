import csv
import os

RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "unified")

# Reference dictionaries

MS_MEDIUM_MG_PER_L = {
    # Macronutrients
    "KNO3": 1900,                # Potassium nitrate
    "NH4NO3": 1650,              # Ammonium nitrate
    "CaCl2.2H2O": 440,           # Calcium chloride dihydrate
    "MgSO4.7H2O": 370,           # Magnesium sulfate heptahydrate
    "KH2PO4": 170,               # Potassium phosphate monobasic
    # Micronutrients
    "MnSO4.4H2O": 22.3,          # Manganese sulfate tetrahydrate
    "ZnSO4.7H2O": 8.6,           # Zinc sulfate heptahydrate
    "H3BO3": 6.2,                # Boric acid
    "KI": 0.83,                  # Potassium iodide
    "CuSO4.5H2O": 0.025,         # Copper(II) sulfate pentahydrate
    "Na2MoO4.2H2O": 0.25,        # Sodium molybdate dihydrate
    "CoCl2.6H2O": 0.025,         # Cobalt chloride hexahydrate
    # Iron source
    "FeSO4.7H2O": 27.8,          # Ferrous sulfate heptahydrate
    "Na2EDTA.2H2O": 37.3,        # EDTA disodium salt, dihydrate
    # Vitamins
    "C6H12O6_myo_inositol": 100, # Myo-inositol
    "C6H5NO2": 0.5,              # Nicotinic acid
    "C8H11NO3.HCl": 0.5,         # Pyridoxine HCl
    "C12H17ClN4OS.HCl": 0.1,     # Thiamine HCl
    "C2H5NO2": 2.0,              # Glycine
}

# Formula string -> {ion_name: ion count per formula}.
# An empty dict means you can assume the compound cannot be broken down
# any further and should remain as is.
COMPOUND_TO_IONS = {
    # Macronutrients
    "KNO3": {"K+": 1, "NO3-": 1},
    "NH4NO3": {"NH4+": 1, "NO3-": 1},
    "CaCl2.2H2O": {"Ca2+": 1, "Cl-": 2},
    "MgSO4.7H2O": {"Mg2+": 1, "SO4_2-": 1},
    "KH2PO4": {"K+": 1, "HPO4_2-": 1},
    # Micronutrients
    "MnSO4.4H2O": {"Mn2+": 1, "SO4_2-": 1},
    "ZnSO4.7H2O": {"Zn2+": 1, "SO4_2-": 1},
    "H3BO3": {"BO3": 1},
    "KI": {"K+": 1, "I-": 1},
    "CuSO4.5H2O": {"Cu2+": 1, "SO4_2-": 1},
    "Na2MoO4.2H2O": {"Na+": 2, "MoO2_2-": 1},
    "CoCl2.6H2O": {"Co2+": 1, "Cl-": 2},
    # Iron source
    "FeSO4.7H2O": {"Fe2+": 1, "SO4_2-": 1},
    "Na2EDTA.2H2O": {"Na+": 2, "EDTA": 1},
    # Vitamins -- none ionize in this dataset
    "C6H12O6_myo_inositol": {},
    "C6H5NO2": {},
    "C8H11NO3.HCl": {},
    "C12H17ClN4OS.HCl": {}, 
    "C2H5NO2": {},
}

# Formula string -> molecular weight (g/mol).
MOLECULAR_WEIGHTS = {
    # Macronutrients
    "KNO3": 101.10,
    "NH4NO3": 80.04,
    "CaCl2.2H2O": 147.01,
    "MgSO4.7H2O": 246.48,
    "KH2PO4": 136.09,
    # Micronutrients
    "MnSO4.4H2O": 223.06,
    "ZnSO4.7H2O": 287.60,
    "H3BO3": 61.84,
    "KI": 166.00,
    "CuSO4.5H2O": 249.69,
    "Na2MoO4.2H2O": 241.96,
    "CoCl2.6H2O": 237.93,
    # Iron source
    "FeSO4.7H2O": 278.02,
    "Na2EDTA.2H2O": 372.24,
    # Vitamins
    "C6H12O6_myo_inositol": 180.16,
    "C6H5NO2": 123.11,
    "C8H11NO3.HCl": 205.64,
    "C12H17ClN4OS.HCl": 337.30,
    "C2H5NO2": 75.07,
}

# Maps the ion-dictionary keys onto the flat column names used in our
# unified schema (matching Study 2's native ion columns).
ION_TO_COLUMN = {
    "NO3-": "no3_mM",
    "NH4+": "nh4_mM",
    "K+": "k_mM",
    "Ca2+": "ca_mM",
    "Mg2+": "mg_mM",
    "SO4_2-": "so4_mM",
    "HPO4_2-": "po4_mM",
    "Cl-": "cl_mM",
}


def compute_ms_ion_concentrations():
    """
    Converts the MS basal medium (given in mg/L) into per-ion mM
    concentrations, via molecular weight -> molarity -> ion dissociation.
    """
    ion_totals_mM = {}
    for compound, mg_per_L in MS_MEDIUM_MG_PER_L.items():
        molw = MOLECULAR_WEIGHTS.get(compound)
        mM = (mg_per_L / molw) 
        for ion, count in COMPOUND_TO_IONS.get(compound, {}).items():
            ion_totals_mM[ion] = ion_totals_mM.get(ion, 0) + mM * count

    result = {}
    for ion, column in ION_TO_COLUMN.items():
        result[column] = round(ion_totals_mM.get(ion, 0), 2)
    return result


MS_FIXED_MACRONUTRIENTS = compute_ms_ion_concentrations()

# Study 2's fixed hormone background\
STUDY2_FIXED_HORMONES_mgL = {"bap_mgL": 1.0, "iba_mgL": 0.1}

# Unified schema
FIELDNAMES = [
    "source_table", "media_id", "varied_input_type",
    "bap_mgL", "kin_mgL", "tdz_mgL", "iba_mgL", "naa_mgL",
    "no3_mM", "nh4_mM", "k_mM", "ca_mM", "mg_mM", "so4_mM", "po4_mM", "cl_mM",
    "ns_mean", "ns_se", "ls_mean", "ls_se",
    "cw_mean", "cw_se", "qi_mean", "qi_se",
]


BLANK_HORMONES = {"bap_mgL": 0, "kin_mgL": 0, "tdz_mgL": 0, "iba_mgL": 0, "naa_mgL": 0}

def load_study1_hormone_table(path, source_table, auxin_column):
    """
    Study 1 tables vary BAP against one auxin (IBA or NAA). Macronutrients
    are fixed at standard MS, filled in with the computed constant.
    """
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f, skipinitialspace=True):
            row = {fn: "" for fn in FIELDNAMES}
            row["source_table"] = source_table
            row["media_id"] = r["media"]
            row["varied_input_type"] = "hormone"

            hormones = BLANK_HORMONES.copy()
            hormones["bap_mgL"] = float(r["bap_mgL"])
            hormones[auxin_column] = float(r[auxin_column])
            row.update(hormones)

            row.update(MS_FIXED_MACRONUTRIENTS)  # fixed background
            for result in ("ns", "ls", "cw", "qi"):
                row[f"{result}_mean"] = float(r[f"{result}_mean"])
                row[f"{result}_se"] = float(r[f"{result}_se"])
            rows.append(row)
    return rows


def load_study2_macronutrient_table(path, source_table):
    """
    Study 2's table varies macronutrient ion concentrations in
    mM. Hormones are fixed at 1 mgL BAP + 0.1 mgL IBA, filled in with
    that constant.
    """
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f, skipinitialspace=True):
            row = {fn: "" for fn in FIELDNAMES}
            row["source_table"] = source_table
            row["media_id"] = r["media"]
            row["varied_input_type"] = "macronutrient"

            hormones = BLANK_HORMONES.copy()
            hormones.update(STUDY2_FIXED_HORMONES_mgL)
            row.update(hormones)

            for ion_col in ("no3_mM", "nh4_mM", "k_mM", "ca_mM", "mg_mM", "so4_mM", "po4_mM", "cl_mM"):
                row[ion_col] = float(r[ion_col])

            for result in ("ns", "ls", "cw", "qi"):
                row[f"{result}_mean"] = float(r[f"{result}_mean"])
                row[f"{result}_se"] = float(r[f"{result}_se"])
            rows.append(row)
    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_rows = []
    all_rows += load_study1_hormone_table(os.path.join(RAW_DIR, "study1_table1.csv"), source_table="study1_table1", auxin_column="iba_mgL")
    all_rows += load_study1_hormone_table(os.path.join(RAW_DIR, "study1_table2.csv"), source_table="study1_table2", auxin_column="naa_mgL")
    all_rows += load_study2_macronutrient_table(os.path.join(RAW_DIR, "study2_table1.csv"), source_table="study2_table1")

    out_path = os.path.join(OUT_DIR, "protocols.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"Wrote {len(all_rows)} rows to {out_path}")

    # --- Sanity check ---------------------------------------------------
    # Study 2's own Table 1 reports a "MS" control row with ion
    # concentrations for standard MS medium. Compare that to what we just
    # computed from MS_MEDIUM_MG_PER_L, as a check on the conversion.
    print("\nSanity check -- computed MS macronutrients (mM) vs. Study 2's own reported MS control row:")
    reported_ms_control = {
        "no3_mM": 39.41, "nh4_mM": 20.62, "k_mM": 20.04, "ca_mM": 3.96,
        "mg_mM": 1.5, "so4_mM": 1.74, "po4_mM": 1.25, "cl_mM": 2.99,
    }
    for col in ION_TO_COLUMN.values():
        computed = MS_FIXED_MACRONUTRIENTS[col]
        reported = reported_ms_control[col]
        flag = "  <-- discrepancy" if abs(computed - reported) > 0.15 else ""
        print(f"  {col}: computed={computed:.2f}  reported={reported:.2f}{flag}")


if __name__ == "__main__":
    main()