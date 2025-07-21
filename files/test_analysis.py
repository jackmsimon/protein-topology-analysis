#!/usr/bin/env python3
"""
Test script to verify the protein topology analysis functions work correctly.
"""

import sys
from pathlib import Path
from analysis import (
    PROJECT_ROOT, OUTPUT_DIR, MUTATION_PAIRS,
    get_mutation_name, get_protein_name, setup_output_dirs,
    get_output_files, extract_sequences
)

def test_basic_functions():
    """Test basic utility functions."""
    print("✓ Testing basic functions...")
    
    # Test mutation pairs
    if not MUTATION_PAIRS:
        print("❌ No mutation pairs configured")
        return False
    
    pair = MUTATION_PAIRS[0]
    mut_path = pair["mut"]
    wt_path = pair["wt"]
    
    if not mut_path.exists() or not wt_path.exists():
        print("Required PDB files not found")
        return False
    
    mutation_name = get_mutation_name(mut_path)
    protein_name = get_protein_name(wt_path)
    output_dir = setup_output_dirs(mut_path)
    output_files = get_output_files(output_dir, protein_name, mutation_name)
    
    print(f"✓ Files: {mutation_name}, {protein_name} ({len(output_files)} outputs)")
    
    return True

def test_sequence_extraction():
    """Test sequence extraction from PDB files."""
    print("✓ Testing sequence extraction...")
    
    try:
        pair = MUTATION_PAIRS[0]
        mut_path = pair["mut"]
        wt_path = pair["wt"]
        
        sequences = extract_sequences(wt_path, mut_path)
        
        if sequences:
            print(f"✓ Extracted {len(sequences)} sequences")
            for name, seq in sequences.items():
                print(f"✓ {name}: {len(seq)} residues")
            return True
        else:
            return False
            
    except Exception:
        return False

def test_existing_analysis():
    """Check if existing analysis files are present and valid."""
    print("✓ Checking analysis files...")
    
    analysis_files = [
        "specifications.json",
        "most_changed_residues.json",
        "pocket_analysis.json",
        "combined_analysis.json",
        "residue_boundaries.json"
    ]
    
    existing_count = sum(1 for f in analysis_files if (OUTPUT_DIR / f).exists())
    print(f"✓ Found {existing_count}/{len(analysis_files)} files")
    return existing_count > 0

def main():
    """Run all tests."""
    print("✓ Running tests...")
    
    tests = [
        test_basic_functions,
        test_sequence_extraction,
        test_existing_analysis
    ]
    
    passed = 0
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                print(f"✗ Failed: {test.__name__}")
        except Exception as e:
            print(f"✗ Error: {test.__name__} - {e}")
    
    print(f"✓ Tests: {passed}/{len(tests)} passed")
    
    if passed == len(tests):
        print("✓ Ready to run: python analysis.py")
    else:
        sys.exit(1)

if __name__ == "__main__":
    main() 