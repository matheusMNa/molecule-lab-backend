from molecule_lab.chem.rupture_database import (
    RuptureDatabase,
)

def test_database_loads():
    db = RuptureDatabase.from_json(
        "tests/data/test_rupture_db.json"
	)
    
    assert db is not None
    
def test_energy_lookup():
    db = RuptureDatabase.from_json(
        "tests/data/test_rupture_db.json"
	)
    
    energy = db.get_energy("CCO")
    
    assert energy == 410.0