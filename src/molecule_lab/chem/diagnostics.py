'''
Mas que química o simulador pensa que existe??
'''

from __future__ import annotations

from typing import Any

def summarize_bond_env(topology) -> list[dict[str, Any]]:
	summary = []
    
	for idx, bond in enumerate(topology.bonds):
		i, j, r0, de, alpha, signature, source = bond
        
		symbols = f'{topology.symbols[i]} - {topology.symbols[j]}'
        
		summary.append(
			{
				"bond_index": idx,
				"atoms": (i, j),
                "bond": symbols,
                "signature": signature,
                "r0": r0,
                "De": de,
                "alpha": alpha,
				"source": source,
			}
		)
	
	return summary