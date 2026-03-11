"""
Thermodynamic On-Target Probe Scoring

"""

import numpy as np
import math
from typing import Dict, List, Optional
from collections import Counter
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
import primer3


class ThermodynamicProbeScorer:
    """
    Thermodynamic-based probe quality scorer.
    Evaluates probes based on intrinsic sequence properties.
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
        Initialize scorer with experimental conditions and design criteria.
        
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
    
    def calculate_on_target_score(self, sequence: str) -> Dict:
        """
        Calculate comprehensive on-target quality score (0-100).
        """
        sequence = sequence.upper()
        
        
        passes, rejection_reason = self.passes_hard_filters(sequence)
        
        if not passes:
            
            return {
                'on_target_score': 0.0,
                'tm': 0.0,
                'gc_content': 0.0,
                'probe_length': len(sequence),
                'has_homopolymer': True,
                'complexity': 0.0,
                'secondary_structure_penalty': 0.0,
                'rejected': True,
                'rejection_reason': rejection_reason,
                'component_scores': {
                    'tm_score': 0.0,
                    'gc_score': 0.0,
                    'length_score': 0.0,
                    'homopolymer_score': 0.0,
                    'secondary_structure_score': 0.0,
                    'complexity_score': 0.0
                }
            }
        
       
        tm = self.calculate_tm(sequence)
        gc_content = self.calculate_gc_content(sequence)
        probe_length = len(sequence)
        has_homopolymer = self.check_homopolymer_runs(sequence)
        complexity = self.calculate_complexity(sequence)
        sec_struct_penalty = self.calculate_secondary_structure_penalty(sequence)
        
        tm_diff = abs(tm - self.target_tm)
        tm_score = max(0, 100 - (tm_diff * 10))  
    
        if self.optimal_gc_range[0] <= gc_content <= self.optimal_gc_range[1]:
            gc_score = 100.0
        else:
            if gc_content < self.optimal_gc_range[0]:
                gc_diff = self.optimal_gc_range[0] - gc_content
            else:
                gc_diff = gc_content - self.optimal_gc_range[1]
            gc_score = max(0, 100 - (gc_diff * 2))
        
        if self.optimal_length_range[0] <= probe_length <= self.optimal_length_range[1]:
            length_score = 100.0
        else:
            if probe_length < self.optimal_length_range[0]:
                length_diff = self.optimal_length_range[0] - probe_length
            else:
                length_diff = probe_length - self.optimal_length_range[1]
            length_score = max(0, 100 - (length_diff * 5))
      
        homopolymer_score = 0.0 if has_homopolymer else 100.0
        
        sec_struct_score = (1.0 - sec_struct_penalty) * 100
        
        complexity_score = complexity * 100
        
        weights = {
            'tm': 0.30,
            'gc': 0.20,
            'length': 0.15,
            'homopolymer': 0.15,
            'secondary_structure': 0.10,
            'complexity': 0.10
        }
      
        on_target_score = (
            weights['tm'] * tm_score +
            weights['gc'] * gc_score +
            weights['length'] * length_score +
            weights['homopolymer'] * homopolymer_score +
            weights['secondary_structure'] * sec_struct_score +
            weights['complexity'] * complexity_score
        )
        
        return {
            'on_target_score': round(on_target_score, 2),
            'tm': round(tm, 2),
            'gc_content': round(gc_content, 2),
            'probe_length': probe_length,
            'has_homopolymer': has_homopolymer,
            'complexity': round(complexity, 3),
            'secondary_structure_penalty': round(sec_struct_penalty, 3),
            'rejected': False,
            'rejection_reason': None,
            'component_scores': {
                'tm_score': round(tm_score, 2),
                'gc_score': round(gc_score, 2),
                'length_score': round(length_score, 2),
                'homopolymer_score': round(homopolymer_score, 2),
                'secondary_structure_score': round(sec_struct_score, 2),
                'complexity_score': round(complexity_score, 2)
            }
        }
    
    def score_probe_set(self, sequences: List[str]) -> List[Dict]:
        """
        Score multiple probes and return sorted results.
        """
        results = []
        
        for seq in sequences:
            score_info = self.calculate_on_target_score(seq)
            score_info['sequence'] = seq
            results.append(score_info)
   
        results.sort(key=lambda x: x['on_target_score'], reverse=True)
        
        return results
    
    def _reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence."""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join(complement.get(base, 'N') for base in reversed(sequence))
