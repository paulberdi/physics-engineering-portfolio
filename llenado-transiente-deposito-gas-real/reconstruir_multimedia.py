#!/usr/bin/env python3
"""Reconstruye los tres archivos multimedia grandes a partir de sus fragmentos.
Uso: python reconstruir_multimedia.py
"""
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
manifests = sorted((BASE / "Proyecto_Final_Corregido").rglob("*.parts/manifest.json"))
if len(manifests) != 3:
    raise SystemExit(f"Se esperaban tres archivos fragmentados; encontrados: {len(manifests)}")

for manifest in manifests:
    meta = json.loads(manifest.read_text(encoding="utf-8"))
    target = manifest.parent.parent / meta["filename"]
    expected_size = meta["size_bytes"]
    expected_sha = meta["sha256"]
    if target.exists():
        if target.stat().st_size == expected_size and hashlib.sha256(target.read_bytes()).hexdigest() == expected_sha:
            print(f"Ya existe y coincide: {target.relative_to(BASE)}")
            continue
        raise SystemExit(f"No se reemplaza un archivo diferente: {target}")
    temp = target.with_suffix(target.suffix + ".partial")
    digest = hashlib.sha256()
    count = 0
    try:
        with temp.open("wb") as output:
            for i in range(1, meta["parts"] + 1):
                part = manifest.parent / f"part{i:03d}"
                data = part.read_bytes()
                output.write(data)
                digest.update(data)
                count += len(data)
        if count != expected_size or digest.hexdigest() != expected_sha:
            raise ValueError(f"Fallo de integridad para {target.name}")
        temp.replace(target)
        print(f"Restaurado y verificado: {target.relative_to(BASE)}")
    finally:
        if temp.exists():
            temp.unlink()
print("Reconstrucción completa.")
