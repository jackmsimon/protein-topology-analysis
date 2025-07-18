# Protein Topology Analysis Pipeline

An automated computational pipeline for analyzing protein mutations through contact map analysis and binding pocket identification. The system takes wildtype and mutant protein structures as input and produces comprehensive structural analyses designed to identify key residues for protein design and drug targeting applications.

The pipeline integrates established computational biology tools including MDTraj for molecular dynamics trajectory analysis, fpocket for binding cavity detection, and advanced contact map algorithms for structural comparison. Large language models then interpret structural data and make intelligent binding site selections, enabling the workflow to operate end-to-end with minimal manual intervention.

This analysis framework addresses a critical need in structural biology and drug discovery, where understanding mutation-induced structural changes is essential for rational protein design. The system identifies residues with significant contact changes between wildtype and mutant structures, maps these changes to druggable binding pockets, and generates prioritized hotspot specifications for downstream protein engineering workflows.

The complete pipeline from PDB input to final hotspot specifications typically completes within 1-2 minutes. 
