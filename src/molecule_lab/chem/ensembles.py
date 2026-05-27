"""
Extração de ambientes locais ao redor de ligações e geração de
assinaturas canônicas para consulta ao banco de rupturas.

5 átomos como limite superior, p/ moleculas menores BFS esgota
todos os átomos disponíveis e retorna o ambiente completo 

 mecanismo de lookup é uniforme.
"""
from __future__ import annotations

from collections import deque

from rdkit import Chem


def extract_local_environment(
        mol: Chem.Mol,
        bond: Chem.Bond,
        size: int = 5,
) -> list[int]:
    """Extrai um ambiente determinístico local ao redor de uma ligação.
    O ambiente sempre contém os dois átomos da ligação;
    expande via busca em largura;
    usa ordenação por índice atômico para garantir determinismo;
    para ao atingir ``size`` átomos ou esgotar a molécula

    Parâmetros
    ----------
    mol:
        Molécula RDKit (com hidrogênios explícitos, se necessário).
    bond:
        Ligação-alvo cujo ambiente será extraído.
    size:
        Número máximo de átomos no ambiente. Padrão: 5.

    Retorno
    -------
    list[int]
        Índices atômicos ordenados que pertencem ao ambiente.
    """
    begin = bond.GetBeginAtomIdx()
    end = bond.GetEndAtomIdx()

    selected = {begin, end}
    queue = deque([begin, end])

    while queue and len(selected) < size:
        current = queue.popleft()
        atom = mol.GetAtomWithIdx(current)

        neighbors = sorted(
            nbr.GetIdx()
            for nbr in atom.GetNeighbors()
        )

        for nbr_idx in neighbors:
            if nbr_idx not in selected:
                selected.add(nbr_idx)
                queue.append(nbr_idx)

            if len(selected) >= size:
                break

    return sorted(selected)


def canonical_environment_signature(
        mol: Chem.Mol,
        atom_indices: list[int],
) -> str:
    """Gera uma assinatura SMILES canônica para um ambiente local.

    A assinatura é usada como chave de lookup no banco de rupturas.
    Ambientes quimicamente equivalentes == msm assinatura

    Parâmetros
    mol:
        Molécula RDKit da qual o fragmento será extraído.
    atom_indices:
        Índices dos átomos que compõem o ambiente.

    Retorno
    str
        Fragmento SMILES canônico do ambiente.
    """
    return Chem.MolFragmentToSmiles(
        mol,
        atomsToUse=sorted(atom_indices),
        canonical=True,
    )