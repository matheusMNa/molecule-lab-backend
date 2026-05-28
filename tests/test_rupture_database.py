from molecule_lab.chem.rupture_database import (
    RuptureDatabase,
)
import pytest

def test_database_loads():
    db = RuptureDatabase.from_json(
        "tests/data/rupture_db.json"
	)
    
    assert db is not None
    
def test_energy_lookup():
    db = RuptureDatabase.from_json(
        "tests/data/rupture_db.json"
	)
    
    energy = db.get_energy("CCO")
    
    assert energy == pytest.approx(88866.2, rel=1e-3)