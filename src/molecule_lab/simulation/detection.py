"""Bond breaking detection via temperature threshold."""

from __future__ import annotations
import math
import numpy as np


def detect_broken_bonds(
    pos: np.ndarray,
    bonds: list[tuple[int, int, float, float, float]],
    current_temperature: float,
    break_distance_factor: float = 1.6,
) -> list[dict]:
    """Detect bonds that have exceeded their rupture temperature threshold.

    The rupture criterion is purely thermal: a bond breaks when the current
    simulation temperature meets or exceeds the threshold temperature stored
    in the bond's ``de`` field (which now holds ``temperatura_K`` from the
    professor's database).

    Bonds with a negative threshold (spontaneously unstable) are always
    considered broken, regardless of the current temperature.

    A secondary distance guard (``break_distance_factor``) is kept to avoid
    flagging bonds that are thermally above threshold but geometrically intact
    — this prevents false positives during the early heating ramp.

    Parameters
    ----------
    pos:
        Atom positions array of shape ``(n_atoms, 3)``.
    bonds:
        List of bond tuples ``(i, j, r0, rupture_temp_K, alpha, ...)``.
        ``rupture_temp_K`` occupies the slot formerly used by Morse ``De``.
    current_temperature:
        Instantaneous kinetic temperature of the simulation in Kelvin.
    break_distance_factor:
        Bond is only declared broken if ``r / r0 >= break_distance_factor``
        *or* the rupture threshold is negative (spontaneous instability).

    Returns
    -------
    list[dict]
        One entry per broken bond with diagnostic fields.
    """
    broken: list[dict] = []

    for idx, (i, j, r0, rupture_temp, alpha, *_) in enumerate(bonds):
        spontaneous = rupture_temp < 0.0
        temp_trigger = spontaneous or (current_temperature >= rupture_temp)

        if not temp_trigger:
            continue

        r = float(np.linalg.norm(pos[i] - pos[j]))
        distance_ratio = r / r0 if r0 > 0 else float("inf")
        distance_trigger = distance_ratio >= break_distance_factor

        # Spontaneously unstable bonds bypass the distance guard.
        if not (distance_trigger or spontaneous):
            continue

        reason = "spontaneous" if spontaneous else "temperature_threshold"

        broken.append(
            {
                "bond_index": idx,
                "i": i,
                "j": j,
                "distance": r,
                "r0": float(r0),
                "distance_ratio": float(distance_ratio),
                "rupture_temp_K": float(rupture_temp),
                "current_temp_K": float(current_temperature),
                # Keep legacy key names so engine._format_broken_bonds works
                # without changes (V and De are repurposed for display).
                "V": float(current_temperature),
                "De": float(rupture_temp),
                "reason": reason,
            }
        )

    return broken
