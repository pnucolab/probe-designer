"""
Single source of truth for SHARP-FISH organism configuration.

Loads `config/organisms.yml` and exposes lookup helpers so backend and
pipeline code never need to hard-code `'human'` / `'mouse'` or genome paths.
Add a new organism by editing the YAML — no code changes required.
"""

from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import yaml


_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "config",
    "organisms.yml",
)

_LOCK = threading.Lock()
_CACHE: Optional["Registry"] = None


@dataclass(frozen=True)
class Microbiome:
    id: str
    display_name: str
    data_dir: str
    source_url: Optional[str] = None


@dataclass(frozen=True)
class Host:
    id: str
    display_name: str
    latin_binomial: Optional[str]
    ncbi_tax_id: Optional[int]
    transcript_id_prefixes: tuple
    genome_dir: Optional[str]
    genome_file_regex: Optional[str]
    transcript_chunks_dir: Optional[str]
    reference_url: Optional[str]
    default: bool = False
    microbiomes: tuple = ()

    @property
    def has_genome(self) -> bool:
        return bool(self.genome_dir and self.genome_file_regex)


@dataclass(frozen=True)
class Registry:
    hosts: tuple  # tuple[Host, ...]
    _hosts_by_id: Dict[str, Host] = field(default_factory=dict, repr=False)
    _microbiomes_by_id: Dict[str, tuple] = field(default_factory=dict, repr=False)
    _host_by_binomial: Dict[str, Host] = field(default_factory=dict, repr=False)
    _host_by_tax_id: Dict[int, Host] = field(default_factory=dict, repr=False)
    _host_by_prefix: Dict[str, Host] = field(default_factory=dict, repr=False)
    _microbiome_to_host: Dict[str, Host] = field(default_factory=dict, repr=False)

    # ---- top-level lookups -------------------------------------------------
    def host_ids(self) -> List[str]:
        return [h.id for h in self.hosts]

    def microbiome_ids(self) -> List[str]:
        return list(self._microbiomes_by_id.keys())

    def get_host(self, host_id: str) -> Optional[Host]:
        return self._hosts_by_id.get(host_id)

    def get_microbiome(self, microbiome_id: str):
        return self._microbiomes_by_id.get(microbiome_id)  # (host, microbiome)

    def host_for_microbiome(self, microbiome_id: str) -> Optional[Host]:
        return self._microbiome_to_host.get(microbiome_id)

    def host_for_binomial(self, binomial: str) -> Optional[Host]:
        if not binomial:
            return None
        return self._host_by_binomial.get(binomial.strip().lower())

    def host_for_tax_id(self, tax_id) -> Optional[Host]:
        try:
            return self._host_by_tax_id.get(int(tax_id))
        except (TypeError, ValueError):
            return None

    def host_for_transcript_id(self, rname: str) -> Optional[Host]:
        if not rname:
            return None
        # Strip Ensembl version suffix
        base = rname.split('.', 1)[0]
        # Iterate longest prefix first so specific prefixes (ENSGALT, 7 chars)
        # match before generic ones (ENSG, 4 chars). Without this, a chicken
        # `ENSGALT…` ID would incorrectly resolve to human via `ENSG`.
        for prefix in sorted(self._host_by_prefix.keys(), key=len, reverse=True):
            if base.startswith(prefix):
                return self._host_by_prefix[prefix]
        return None

    def is_host_species(self, species_id: str) -> bool:
        return species_id in self._hosts_by_id

    def is_microbiome_species(self, species_id: str) -> bool:
        return species_id in self._microbiomes_by_id

    def default_host(self) -> Host:
        for h in self.hosts:
            if h.default:
                return h
        return self.hosts[0]

    def to_json_dict(self) -> dict:
        """Shape suitable for serving to the frontend via JSON."""
        return {
            "hosts": [
                {
                    "id": h.id,
                    "display_name": h.display_name,
                    "latin_binomial": h.latin_binomial,
                    "default": h.default,
                    "has_genome": h.has_genome,
                    "reference_url": h.reference_url,
                    "transcript_id_prefixes": list(h.transcript_id_prefixes),
                    "microbiomes": [
                        {
                            "id": m.id,
                            "display_name": m.display_name,
                            "source_url": m.source_url,
                        }
                        for m in h.microbiomes
                    ],
                }
                for h in self.hosts
            ],
        }


def _build_registry(data: dict) -> Registry:
    hosts: List[Host] = []
    hosts_by_id: Dict[str, Host] = {}
    microbiomes_by_id: Dict[str, tuple] = {}
    host_by_binomial: Dict[str, Host] = {}
    host_by_tax_id: Dict[int, Host] = {}
    host_by_prefix: Dict[str, Host] = {}
    microbiome_to_host: Dict[str, Host] = {}

    for raw in data.get("hosts", []):
        microbiomes = tuple(
            Microbiome(
                id=m["id"],
                display_name=m["display_name"],
                data_dir=m["data_dir"],
                source_url=m.get("source_url"),
            )
            for m in raw.get("microbiomes", []) or []
        )
        host = Host(
            id=raw["id"],
            display_name=raw["display_name"],
            latin_binomial=raw.get("latin_binomial"),
            ncbi_tax_id=raw.get("ncbi_tax_id"),
            transcript_id_prefixes=tuple(raw.get("transcript_id_prefixes") or []),
            genome_dir=raw.get("genome_dir"),
            genome_file_regex=raw.get("genome_file_regex"),
            transcript_chunks_dir=raw.get("transcript_chunks_dir"),
            reference_url=raw.get("reference_url"),
            default=bool(raw.get("default", False)),
            microbiomes=microbiomes,
        )
        hosts.append(host)
        hosts_by_id[host.id] = host
        if host.latin_binomial:
            host_by_binomial[host.latin_binomial.strip().lower()] = host
        if host.ncbi_tax_id is not None:
            host_by_tax_id[int(host.ncbi_tax_id)] = host
        for prefix in host.transcript_id_prefixes:
            host_by_prefix[prefix] = host
        for m in microbiomes:
            microbiomes_by_id[m.id] = (host, m)
            microbiome_to_host[m.id] = host

    return Registry(
        hosts=tuple(hosts),
        _hosts_by_id=hosts_by_id,
        _microbiomes_by_id=microbiomes_by_id,
        _host_by_binomial=host_by_binomial,
        _host_by_tax_id=host_by_tax_id,
        _host_by_prefix=host_by_prefix,
        _microbiome_to_host=microbiome_to_host,
    )


def load_registry(path: Optional[str] = None, force_reload: bool = False) -> Registry:
    """Return the cached registry, parsing the YAML on first access."""
    global _CACHE
    with _LOCK:
        if _CACHE is not None and not force_reload:
            return _CACHE
        cfg_path = path or _CONFIG_PATH
        with open(cfg_path, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
        _CACHE = _build_registry(data)
        return _CACHE


# Convenience top-level wrappers --------------------------------------------

def get_host(host_id: str):
    return load_registry().get_host(host_id)


def get_microbiome(microbiome_id: str):
    return load_registry().get_microbiome(microbiome_id)


def host_for_microbiome(microbiome_id: str):
    return load_registry().host_for_microbiome(microbiome_id)


def host_for_binomial(binomial: str):
    return load_registry().host_for_binomial(binomial)


def host_for_tax_id(tax_id):
    return load_registry().host_for_tax_id(tax_id)


def host_for_transcript_id(rname: str):
    return load_registry().host_for_transcript_id(rname)


def is_host_species(species_id: str) -> bool:
    return load_registry().is_host_species(species_id)


def is_microbiome_species(species_id: str) -> bool:
    return load_registry().is_microbiome_species(species_id)


def host_ids():
    return load_registry().host_ids()


def microbiome_ids():
    return load_registry().microbiome_ids()


def default_host():
    return load_registry().default_host()
