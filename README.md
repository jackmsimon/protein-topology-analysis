# Protein Topology Analysis Pipeline

An automated computational pipeline for analyzing protein mutations through contact map analysis and binding pocket identification. The system takes wildtype and mutant protein structures as input and produces comprehensive structural analyses designed to identify key residues for protein design and drug targeting applications.

The pipeline integrates established computational biology tools including MDTraj for molecular dynamics trajectory analysis, fpocket for binding cavity detection, and contact mapping algorithms for structural comparison. Large language models then interpret structural data and make intelligent binding site selections, enabling the workflow to operate end-to-end with minimal manual intervention.

The complete pipeline from PDB input to final hotspot specifications typically completes within 1-2 minutes. 
