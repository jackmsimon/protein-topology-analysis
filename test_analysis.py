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
    print("=== Testing Basic Functions ===")
    
    # Test path setup
    print(f"✓ PROJECT_ROOT: {PROJECT_ROOT}")
    print(f"✓ OUTPUT_DIR: {OUTPUT_DIR}")
    
    # Test mutation pairs
    if not MUTATION_PAIRS:
        print("❌ No mutation pairs configured")
        return False
    
    pair = MUTATION_PAIRS[0]
    mut_path = pair["mut"]
    wt_path = pair["wt"]
    
    print(f"✓ Mutation file: {mut_path} (exists: {mut_path.exists()})")
    print(f"✓ Wildtype file: {wt_path} (exists: {wt_path.exists()})")
    
    if not mut_path.exists() or not wt_path.exists():
        print("❌ Required PDB files not found")
        return False
    
    # Test name extraction
    mutation_name = get_mutation_name(mut_path)
    protein_name = get_protein_name(wt_path)
    print(f"✓ Mutation name: {mutation_name}")
    print(f"✓ Protein name: {protein_name}")
    
    # Test directory setup
    output_dir = setup_output_dirs(mut_path)
    print(f"✓ Output directory: {output_dir}")
    
    # Test output file paths
    output_files = get_output_files(output_dir, protein_name, mutation_name)
    print(f"✓ Output files configured: {len(output_files)} files")
    
    return True

def test_sequence_extraction():
    """Test sequence extraction from PDB files."""
    print("\n=== Testing Sequence Extraction ===")
    
    try:
        pair = MUTATION_PAIRS[0]
        mut_path = pair["mut"]
        wt_path = pair["wt"]
        
        sequences = extract_sequences(wt_path, mut_path)
        
        if sequences:
            print(f"✓ Extracted {len(sequences)} sequences:")
            for name, seq in sequences.items():
                print(f"  - {name}: {len(seq)} residues")
                if len(seq) > 50:
                    print(f"    Preview: {seq[:50]}...")
                else:
                    print(f"    Sequence: {seq}")
            return True
        else:
            print("❌ No sequences extracted")
            return False
            
    except Exception as e:
        print(f"❌ Error in sequence extraction: {e}")
        return False

def test_existing_analysis():
    """Check if existing analysis files are present and valid."""
    print("\n=== Checking Existing Analysis Files ===")
    
    analysis_files = [
        "specifications.json",
        "most_changed_residues.json",
        "pocket_analysis.json",
        "combined_analysis.json",
        "residue_boundaries.json"
    ]
    
    existing_count = 0
    for filename in analysis_files:
        filepath = OUTPUT_DIR / filename
        if filepath.exists():
            print(f"✓ {filename} exists ({filepath.stat().st_size} bytes)")
            existing_count += 1
        else:
            print(f"- {filename} not found")
    
    print(f"✓ Found {existing_count}/{len(analysis_files)} analysis files")
    return existing_count > 0

def main():
    """Run all tests."""
    print("Protein Topology Analysis - Test Suite")
    print("=" * 50)
    
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
                print(f"❌ {test.__name__} failed")
        except Exception as e:
            print(f"❌ {test.__name__} error: {e}")
    
    print(f"\n=== Test Results ===")
    print(f"Passed: {passed}/{len(tests)} tests")
    
    if passed == len(tests):
        print("✓ All tests passed! The repository is ready to use.")
        print("\nTo run the full analysis:")
        print("  python analysis.py")
    else:
        print("❌ Some tests failed. Check the configuration.")
        sys.exit(1)

if __name__ == "__main__":
    main() 