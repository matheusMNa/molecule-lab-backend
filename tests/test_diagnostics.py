from molecule_lab.chem.diagnostics import (
    summarize_bond_env,
)

from molecule_lab.simulation.topology import (
    build_molecule,
    build_topology,
)


def test_bond_environment_summary() -> None:
    mol = build_molecule("C", seed=7)
    topology = build_topology(mol)

    summary = summarize_bond_env(topology)

    assert len(summary) == 4

    first = summary[0]

    assert first["source"] in ("database", "fallback")
    assert "signature" in first
    assert "De" in first
    assert "bond" in first