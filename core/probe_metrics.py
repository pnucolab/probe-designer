"""
Per-probe thermodynamic and sequence metrics.

Measures intrinsic properties of a probe sequence (Tm, GC%, length,
homopolymer runs, Shannon complexity, secondary-structure penalty). No
composite score or ranking is produced — callers apply their own thresholds.
"""

import math
from typing import Dict, List, Optional
from collections import Counter
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import primer3


class ProbeMetricsCalculator:
    """
    Calculates thermodynamic and sequence metrics for individual probes.
    Reports measured properties only; it does not rank or score them.
    """

    def __init__(self,
                 temperature_celsius: float = 37.0,
                 formamide_percent: float = 50.0,
                 na_concentration_mM: float = 390.0,
                 dnac1_nM: float = 25.0,
                 dnac2_nM: float = 25.0,
                 target_tm: float = 47.0,
                 optimal_gc_range: tuple = (40, 60),
                 optimal_length_range: tuple = (30, 37),
                 max_homopolymer: int = 5,
                 enable_hard_filters: bool = True):
        """
        Initialize the calculator with experimental conditions.

        Args:
            temperature_celsius: Hybridization temperature (default 37°C)
            formamide_percent: Formamide concentration (default 50%)
            na_concentration_mM: Sodium concentration in mM (default 390)
            dnac1_nM: Concentration of higher concentration strand in nM (default 25)
            dnac2_nM: Concentration of lower concentration strand in nM (default 25)
            target_tm: Target melting temperature (default 47°C)
            optimal_gc_range: Optimal GC% range (default 40-60%)
            optimal_length_range: Optimal probe length (default 30-37 nt)
            max_homopolymer: Maximum allowed homopolymer run (default 5)
        """
        self.T_celsius = temperature_celsius
        self.formamide = formamide_percent
        self.na_conc_mM = na_concentration_mM
        self.dnac1 = dnac1_nM
        self.dnac2 = dnac2_nM
        self.target_tm = target_tm
        self.optimal_gc_range = optimal_gc_range
        self.optimal_length_range = optimal_length_range
        self.max_homopolymer = max_homopolymer
        self.enable_hard_filters = enable_hard_filters

    def calculate_tm(self, sequence: str) -> float:
        """
        Calculate melting temperature using BioPython.
        """
        tmval = float(('%0.2f' % mt.Tm_NN(sequence, Na=self.na_conc_mM,
                                          dnac1=self.dnac1, dnac2=self.dnac2)))
        tm_corrected = float(('%0.2f' % mt.chem_correction(tmval, fmd=self.formamide)))
        return tm_corrected

    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content percentage using BioPython.
        """
        return gc_fraction(sequence) * 100

    def check_homopolymer_runs(self, sequence: str, max_run: int = None) -> bool:
        """
        Check for problematic homopolymer runs (AAAAA, TTTTT, etc.)

        """
        if max_run is None:
            max_run = self.max_homopolymer

        sequence = sequence.upper()
        for base in ['A', 'T', 'G', 'C']:
            pattern = base * max_run
            if pattern in sequence:
                return True
        return False

    def calculate_complexity(self, sequence: str) -> float:
        """
        Calculate sequence complexity using Shannon entropy.
        """
        sequence = sequence.upper()
        if len(sequence) == 0:
            return 0.0

        counts = Counter(sequence)
        entropy = 0.0

        for base in ['A', 'T', 'G', 'C']:
            if base in counts:
                p = counts[base] / len(sequence)
                if p > 0:
                    entropy -= p * math.log2(p)
        return entropy / 2.0

    def calculate_secondary_structure_penalty(self, sequence: str) -> float:
        """
        Calculate secondary structure penalty using primer3 thermodynamic engine.
        Evaluates hairpin and self-dimer formation potential.
        Returns a penalty between 0.0 (no structure) and 1.0 (strong structure).
        """
        sequence = sequence.upper()

        hairpin = primer3.calc_hairpin(
            sequence,
            mv_conc=self.na_conc_mM,
            dv_conc=0,
            dntp_conc=0,
            dna_conc=self.dnac1,
            temp_c=self.T_celsius
        )

        homodimer = primer3.calc_homodimer(
            sequence,
            mv_conc=self.na_conc_mM,
            dv_conc=0,
            dntp_conc=0,
            dna_conc=self.dnac1,
            temp_c=self.T_celsius
        )

        # delta G in cal/mol, convert to kcal/mol
        hairpin_dg = hairpin.dg / 1000.0
        homodimer_dg = homodimer.dg / 1000.0
        worst_dg = min(hairpin_dg, homodimer_dg)

        # More negative dG = more stable secondary structure = higher penalty
        # Scale: 0 kcal/mol -> 0.0 penalty, -10 kcal/mol or worse -> 1.0 penalty
        penalty = min(max(0, -worst_dg / 10.0), 1.0)
        return penalty

    def passes_hard_filters(self, sequence: str) -> tuple[bool, Optional[str]]:
        """
        Check if sequence passes hard filters (OligoMiner-style automatic rejection).
        """
        if not self.enable_hard_filters:
            return True, None

        if self.check_homopolymer_runs(sequence):
            return False, f"Contains homopolymer run ≥{self.max_homopolymer}bp"


        return True, None

    def calculate_probe_metrics(self, sequence: str) -> Dict:
        """
        Calculate probe quality metrics.
        """
        sequence = sequence.upper()

        passes, rejection_reason = self.passes_hard_filters(sequence)

        if not passes:
            return {
                'tm': 0.0,
                'gc_content': 0.0,
                'probe_length': len(sequence),
                'has_homopolymer': True,
                'complexity': 0.0,
                'secondary_structure_penalty': 0.0,
                'rejected': True,
                'rejection_reason': rejection_reason,
            }

        tm = self.calculate_tm(sequence)
        gc_content = self.calculate_gc_content(sequence)
        probe_length = len(sequence)
        has_homopolymer = self.check_homopolymer_runs(sequence)
        complexity = self.calculate_complexity(sequence)
        sec_struct_penalty = self.calculate_secondary_structure_penalty(sequence)

        return {
            'tm': round(tm, 2),
            'gc_content': round(gc_content, 2),
            'probe_length': probe_length,
            'has_homopolymer': has_homopolymer,
            'complexity': round(complexity, 3),
            'secondary_structure_penalty': round(sec_struct_penalty, 3),
            'rejected': False,
            'rejection_reason': None,
        }

    def calculate_metrics_for_set(self, sequences: List[str]) -> List[Dict]:
        """
        Calculate metrics for multiple probes, preserving input order.
        """
        results = []

        for seq in sequences:
            metrics = self.calculate_probe_metrics(seq)
            metrics['sequence'] = seq
            results.append(metrics)

        return results

    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence))
