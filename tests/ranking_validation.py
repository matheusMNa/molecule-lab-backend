# ============================================================
# SCIENTIFIC RANKING VALIDATION
# ============================================================

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import pandas as pd

# API oficial da simulação
from molecule_lab.simulation import run_simulation


# ============================================================
# CONFIGURAÇÕES
# ============================================================

DATASET_PATH = "data/validation_molecules.csv"

PRESET_NAME = "fast"

SEED = 42

STABLE_TEMPERATURE = 11000.0

OUTPUT_DIR = Path("results")


# ============================================================
# CARREGAR DATASET
# ============================================================

def load_dataset(path: str | Path):

    df = pd.read_csv(path)

    return df.to_dict(orient="records")


# ============================================================
# EXECUTAR SIMULAÇÃO
# ============================================================

def simulate_molecule(
    smiles: str,
    preset_name: str = PRESET_NAME,
    seed: int = SEED,
):

    final_result = None

    for event in run_simulation(
        smiles=smiles,
        preset_name=preset_name,
        seed=seed,
    ):

        if event.event == "result":

            final_result = event.payload

    return final_result


# ============================================================
# EXTRAIR TEMPERATURA
# ============================================================

def extract_break_temperature(
    result_payload: dict,
):

    # molécula quebrou
    if result_payload["result"] == "break":

        return (
            float(
                result_payload[
                    "current_break_temperature"
                ]
            ),
            "break",
        )

    # molécula estável
    return (
        STABLE_TEMPERATURE,
        "stable",
    )


# ============================================================
# VALIDAR MOLÉCULA
# ============================================================

def validate_entry(entry: dict):

    name = entry["name"]

    smiles = entry["smiles"]

    reference_value = float(
        entry["reference_value"]
    )

    print(
        f"Simulando: {name} ({smiles})"
    )

    try:

        result = simulate_molecule(smiles)

        predicted_temp, status = (
            extract_break_temperature(result)
        )

        return {

            "name": name,

            "smiles": smiles,

            "reference_value": reference_value,

            "predicted_temperature": predicted_temp,

            "result": status,
        }

    except Exception as exc:

        return {

            "name": name,

            "smiles": smiles,

            "reference_value": reference_value,

            "predicted_temperature": float("nan"),

            "result": "error",

            "error": str(exc),
        }


# ============================================================
# ADICIONAR RANKS
# ============================================================

def add_ranks(df: pd.DataFrame):

    ranked = df.copy()

    # maior estabilidade = maior temperatura
    ranked["reference_rank"] = (
        ranked["reference_value"]
        .rank(
            ascending=False,
            method="min",
        )
        .astype(int)
    )

    ranked["predicted_rank"] = (
        ranked["predicted_temperature"]
        .rank(
            ascending=False,
            method="min",
        )
        .astype(int)
    )

    ranked["rank_difference"] = (
        ranked["predicted_rank"]
        -
        ranked["reference_rank"]
    )

    return ranked


# ============================================================
# COMPARAÇÃO PAR-A-PAR
# ============================================================

def compare_rankings(df: pd.DataFrame):

    comparable_pairs = 0

    correct_pairs = 0

    wrong_pairs = 0

    rows = df.to_dict(
        orient="records"
    )

    for row_a, row_b in combinations(
        rows,
        2,
    ):

        comparable_pairs += 1

        # ordem experimental
        exp_order = (
            row_a["reference_value"]
            >
            row_b["reference_value"]
        )

        # ordem simulada
        pred_order = (
            row_a["predicted_temperature"]
            >
            row_b["predicted_temperature"]
        )

        if exp_order == pred_order:

            correct_pairs += 1

        else:

            wrong_pairs += 1

    ranking_difference = (
        wrong_pairs /
        comparable_pairs
    ) * 100

    ranking_similarity = (
        correct_pairs /
        comparable_pairs
    ) * 100

    return {

        "n_comparable_pairs": comparable_pairs,

        "n_correct_pairs": correct_pairs,

        "n_wrong_pairs": wrong_pairs,

        "ranking_difference": ranking_difference,

        "ranking_similarity": ranking_similarity,
    }


# ============================================================
# CORRELAÇÕES
# ============================================================

def compute_correlations(
    df: pd.DataFrame,
):

    try:

        from scipy.stats import (
            spearmanr,
            kendalltau,
        )

        spearman = spearmanr(
            df["reference_value"],
            df["predicted_temperature"],
        ).statistic

        kendall = kendalltau(
            df["reference_value"],
            df["predicted_temperature"],
        ).statistic

        return (

            float(spearman),

            float(kendall),
        )

    except Exception:

        return (
            float("nan"),
            float("nan"),
        )


# ============================================================
# EXECUÇÃO PRINCIPAL
# ============================================================

def main():

    print("\n=====================================")
    print("SCIENTIFIC RANKING VALIDATION")
    print("=====================================\n")

    dataset = load_dataset(
        DATASET_PATH
    )

    results = []

    # simulações
    for entry in dataset:

        result = validate_entry(entry)

        results.append(result)

    # dataframe
    df = pd.DataFrame(results)

    # adiciona rankings
    df = add_ranks(df)

    # ordena por estabilidade
    df = df.sort_values(
        by="reference_rank"
    )

    print("\n=====================================")
    print("VALIDATION RESULTS")
    print("=====================================\n")

    display_columns = [

        "name",

        "reference_value",

        "predicted_temperature",

        "reference_rank",

        "predicted_rank",

        "rank_difference",
    ]

    print(
        df[
            display_columns
        ].to_string(index=True)
    )

    # métricas
    summary = compare_rankings(df)

    spearman, kendall = (
        compute_correlations(df)
    )

    print("\n=====================================")
    print("RANKING METRICS")
    print("=====================================\n")

    print(
        f"Pares comparáveis: "
        f"{summary['n_comparable_pairs']}"
    )

    print(
        f"Pares corretos: "
        f"{summary['n_correct_pairs']}"
    )

    print(
        f"Pares invertidos: "
        f"{summary['n_wrong_pairs']}"
    )

    print()

    print(
        f"Ranking difference: "
        f"{summary['ranking_difference']:.2f}%"
    )

    print(
        f"Ranking similarity: "
        f"{summary['ranking_similarity']:.2f}%"
    )

    print()

    print(
        f"Spearman: "
        f"{spearman:.4f}"
    )

    print(
        f"Kendall tau: "
        f"{kendall:.4f}"
    )

    # salvar resultados
    OUTPUT_DIR.mkdir(
        exist_ok=True
    )

    output_csv = (
        OUTPUT_DIR /
        "validation_results.csv"
    )

    df.to_csv(
        output_csv,
        index=False,
    )

    print("\nResultados salvos em:")

    print(
        OUTPUT_DIR.resolve()
    )


# ============================================================
# ENTRYPOINT
# ============================================================

if __name__ == "__main__":

    main()
