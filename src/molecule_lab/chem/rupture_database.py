'''
Esse módulo deve:

carregar as energias de ruptura de ligação
normalizar as assinaturas
prover funções de busca
armazenar informação no cache

provavelmente otimização é armazenar como JSON
'''

from __future__ import annotations

import json
from pathlib import Path

class RuptureDatabase:
    def __init__(self, data: dict[str, float]) -> None:
        self.__data = data
        
    @classmethod
    
    def from_json(cls, path: str | Path) -> "RuptureDatabase":
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return cls(data)
    
    def get_energy(self, signature : str) -> float:
        if signature not in self.__data:
            raise KeyError(
                f"Sem energia de quebra de ligação para a assinatura: {signature}"
			)
        
        return self.__data[signature]