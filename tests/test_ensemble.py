from rdkit import Chem

from molecule_lab.chem.ensembles import (
    extract_local_environment,
    canonical_environment_signature,
)

def test_env_contains_bond_atoms():
    mol = Chem.MolFromSmiles("CCO")
    bond = mol.GetBondWithIdx(1)
    env = extract_local_environment(mol, bond)
    begin = bond.GetBeginAtomIdx()
    end = bond.GetEndAtomIdx() 
    assert begin in env
    assert end in env
    
def test_env_size_limit():
	mol = Chem.MolFromSmiles("CCCCCCCC")
	bond = mol.GetBondWithIdx(2)

	env = extract_local_environment(mol, bond)

	assert len(env) <= 5


def test_env_is_deterministic():
	mol = Chem.MolFromSmiles("CCCO")

	bond = mol.GetBondWithIdx(1)

	env1 = extract_local_environment(mol, bond)
	env2 = extract_local_environment(mol, bond)

	assert env1 == env2

def test_signature_is_deterministic():
	mol = Chem.MolFromSmiles("CCCO")

	bond = mol.GetBondWithIdx(1)

	env = extract_local_environment(mol, bond)

	sig1 =canonical_environment_signature(mol, env)
	sig2 =canonical_environment_signature(mol, env)

	assert sig1 == sig2

def test_equivalent_enviroments_signature():
	
	mol1 = Chem.MolFromSmiles("CCO")
	mol2 = Chem.MolFromSmiles("OCC")
	
	bond1 = mol1.GetBondWithIdx(1)
	bond2 = mol2.GetBondWithIdx(0)
	
	env1 = extract_local_environment(mol1, bond1)
	env2 = extract_local_environment(mol2, bond2)
	
	sig1 = canonical_environment_signature(mol1, env1)
	sig2 = canonical_environment_signature(mol2, env2)
	
	assert sig1 == sig2