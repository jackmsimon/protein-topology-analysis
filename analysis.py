import json
import pandas as pd
import mdtraj as md
import os
import matplotlib.pyplot as plt
from contact_map import ContactFrequency, AtomMismatchedContactDifference
import subprocess
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio.PDB import *
from Bio import PDB
from Bio.SeqUtils import seq1
import anthropic
import warnings
from datetime import datetime

# Suppress matplotlib warnings about contact map visualization
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

# Initialize Anthropic client (you'll need to set your API key)
try:
    client = anthropic.Anthropic(
        api_key="" # Anthropic API key here
    )
except Exception as e:
    client = None

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent

# Define mutation pairs using the current project structure
MUTATION_PAIRS = [
    {
        "mut": PROJECT_ROOT / "files" / "S252W.pdb",
        "wt": PROJECT_ROOT / "files" / "FGFR2.pdb"
    }
]

# Output directory
OUTPUT_DIR = PROJECT_ROOT / "analysis"
OUTPUT_DIR.mkdir(exist_ok=True)

def get_chain_id(chain):
    """Helper function to get chain identifier"""
    if hasattr(chain, 'segment_id') and chain.segment_id:
        return chain.segment_id
    if hasattr(chain, 'chain_id') and chain.chain_id:
        return chain.chain_id
    return str(chain.index)

# Fpocket analysis functions
def parse_pdb_line(line):
    """Parse a PDB file line and return atomic information."""
    if line.startswith(("ATOM", "HETATM")):
        return {
            'record_type': line[0:6].strip(),
            'atom_num': int(line[6:11]),
            'atom_name': line[12:16].strip(),
            'res_name': line[17:20].strip(),
            'chain_id': line[21],
            'res_num': int(line[22:26]),
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54])
        }
    return None

def distance(p1, p2):
    """Calculate Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

def parse_pocket_info(info_text):
    """Parse fpocket info file."""
    pockets = {}
    current_pocket = None
    
    for line in info_text.split('\n'):
        line = line.strip()
        if line.startswith('Pocket'):
            if ':' in line:
                current_pocket = line.split(':')[0].strip()
                pocket_num = int(current_pocket.split()[-1])
                pockets[pocket_num] = {'header': current_pocket, 'info': []}
        elif line and current_pocket and ':' in line:
            key, value = [part.strip().replace('\t', '') for part in line.split(':')]
            cleaned_line = f"{key} : {value}"
            pockets[int(current_pocket.split()[-1])]['info'].append(cleaned_line)
            
    return pockets

def find_pocket_residues(pdb_file, cutoff=4.0):
    """Find residues lining pockets within cutoff distance."""
    pocket_points = defaultdict(list)
    protein_atoms = []
    pocket_residues = defaultdict(set)
    
    with open(pdb_file, 'r') as f:
        for line in f:
            entry = parse_pdb_line(line)
            if not entry:
                continue
                
            if entry['record_type'] == "HETATM":
                if ("POL STP" in line or "APOL STP" in line):
                    pocket_num = entry['res_num']
                    pocket_points[pocket_num].append((entry['x'], entry['y'], entry['z']))
            
            elif entry['record_type'] == "ATOM":
                protein_atoms.append({
                    'res_name': entry['res_name'],
                    'res_num': entry['res_num'],
                    'chain_id': entry['chain_id'],
                    'coords': (entry['x'], entry['y'], entry['z'])
                })

    for pocket_num, points in pocket_points.items():
        for atom in protein_atoms:
            for point in points:
                if distance(atom['coords'], point) <= cutoff:
                    residue = (atom['res_name'], atom['res_num'], atom['chain_id'])
                    pocket_residues[pocket_num].add(residue)
                    break
                
    return {k: sorted(v, key=lambda x: (x[2], x[1])) for k, v in pocket_residues.items()}

def analyze_pockets(pdb_path):
    """Run fpocket analysis on the structure."""
    try:
        base_name = Path(pdb_path).stem
        
        # Create fpocket output directory in the current analysis folder
        fpocket_dir = OUTPUT_DIR / f"{base_name}_out"
        fpocket_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a copy of the PDB file in the analysis directory if needed
        pdb_copy = OUTPUT_DIR / f"{base_name}.pdb"
        if not pdb_copy.exists():
            import shutil
            shutil.copy2(pdb_path, pdb_copy)
        result = subprocess.run(
            ['fpocket', '-f', str(pdb_copy)],
            capture_output=True,
            text=True,
            check=True
        )
        
        # fpocket creates output in the same directory as the input PDB
        info_path = fpocket_dir / f"{base_name}_info.txt"
        pocket_pdb = fpocket_dir / f"{base_name}_out.pdb"
        
        if not info_path.exists():
            raise FileNotFoundError(f"No fpocket output found at {info_path}")
        
        with open(info_path, 'r') as f:
            pocket_info = parse_pocket_info(f.read())
        
        pocket_residues = find_pocket_residues(str(pocket_pdb))
        
        all_pockets = {}
        for pocket_num in sorted(pocket_info.keys()):
            # Use structure name in pocket key format
            pocket_key = f"{base_name}_pocket_{pocket_num}"
            
            formatted_residues = [
                f"{res[0]} {res[2]}{res[1]}"
                for res in pocket_residues.get(pocket_num, [])
            ]
            
            all_pockets[pocket_key] = {
                "info": pocket_info[pocket_num]['info'],
                "residues": sorted(formatted_residues, key=lambda x: int(''.join(filter(str.isdigit, x))))
            }
        
        return all_pockets
        
    except (subprocess.CalledProcessError, Exception):
        return {}

def get_residue_boundaries(pdb_path: str) -> dict:
    """Get first and last residue numbers for chain A only."""
    base_name = Path(pdb_path).stem
    boundaries = {'A': {'first': None, 'last': None}}
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == 'A':  # Only process chain A
                res_num = int(line[22:26])
                
                if boundaries['A']['first'] is None:
                    boundaries['A']['first'] = res_num
                boundaries['A']['last'] = res_num
    
    # Return results with structure name
    return {f"{base_name}": boundaries}

def extract_sequence_from_pdb(pdb_file, chain_id=None):
    """
    Extract amino acid sequence from a PDB file.
    If chain_id is specified, only extract sequence from that chain.
    """
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
        
        # Extract sequence
        sequence = ""
        for model in structure:
            for chain in model:
                # Skip chains that don't match the specified chain_id (if provided)
                if chain_id and chain.id != chain_id:
                    continue
                    
                chain_sequence = ""
                for residue in chain:
                    # Check if it's a standard amino acid residue
                    resname = residue.get_resname()
                    if resname in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
                        try:
                            chain_sequence += seq1(resname)
                        except KeyError:
                            chain_sequence += 'X'  # For non-standard amino acids
                
                # Add chain sequence to total sequence
                sequence += chain_sequence
                
                # If we only want one specific chain, we can break once we've processed it
                if chain_id:
                    break
        
        return sequence
    except Exception:
        return None

def get_mutation_name(mut_path):
    """Extract mutation name from path."""
    return Path(mut_path).stem

def get_protein_name(wt_path):
    """Extract protein name from wildtype path."""
    return Path(wt_path).stem

def setup_output_dirs(mut_path):
    """Setup output directories for a mutation."""
    # Use the global OUTPUT_DIR
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUTPUT_DIR

def get_output_files(output_dir, protein_name, mutation_name):
    """Get output file paths for a mutation."""
    return {
        "diff_csv": output_dir / "residue_contact_difference.csv",
        "output_json": output_dir / "most_changed_residues.json",
        "pocket_json": output_dir / "pocket_analysis.json",
        "boundaries_json": output_dir / "residue_boundaries.json",
        "sequences_json": output_dir / "sequences.json",
        "output_txt": output_dir / "top_contact_diffs.txt",
        "output_png": output_dir / f"{protein_name}_{mutation_name}_contact_diff.png",
        "combined_json": output_dir / "combined_analysis.json",
        "spec_json": output_dir / "specifications.json"
    }

def extract_sequences(wt_path, mut_path):
    """Extract sequences from PDB files and save to JSON."""
    
    # Dictionary to store sequences
    sequences = {}
    
    # Extract sequence from wildtype PDB (all chains)
    wt_seq = extract_sequence_from_pdb(wt_path)
    if wt_seq:
        protein_name = get_protein_name(wt_path)
        sequences[f"{protein_name}_WT"] = wt_seq
    
    # Extract sequence from mutant PDB (all chains)
    mut_seq = extract_sequence_from_pdb(mut_path)
    if mut_seq:
        mutation_name = get_mutation_name(mut_path)
        sequences[f"{protein_name}_{mutation_name}"] = mut_seq
    
    return sequences

def generate_specifications(analysis_data):
    """Generate protein design specifications using Claude AI."""
    if not client:
        print("⚠️ Anthropic client not available, skipping specification generation")
        return {"hotspots": ""}
        
    print("✓ Initializing Anthropic API call...")
    try:
        # Try to load instructions, but use default if not found
        instructions_path = PROJECT_ROOT / "instructions.txt"
        if instructions_path.exists():
            with open(instructions_path, "r") as f:
                instructions = f.read()
            instructions_section = f"""
        Here are instructions for how to interpret and create binding pocket and hotspot specifications:
        --- BEGIN INSTRUCTION ---
        {instructions}
        --- END INSTRUCTION ---"""
        else:
            print("  - No instructions.txt found, using default criteria")
            instructions_section = ""
            
        specifications_prompt = f"""
        Below is a detailed scientific analysis about protein structure and design:
        --- BEGIN ANALYSIS ---
        {json.dumps(analysis_data, indent=2)}
        --- END ANALYSIS ---
        {instructions_section}
        
        Based on the analysis data, identify the most important residues that should be targeted for protein design.
        Output ONLY a JSON object with the following format:
        {{
        "hotspots": "A25,A26,A27,A28,A29,A30,A31,A32,A33,A34"
        }}
        
        CRITICAL SELECTION CRITERIA:
        1. First, identify pockets that contain AT LEAST ONE residue with high contact changes (total_change > 3.0) 
           AND that is also present in interface_residues. These pockets are your primary candidates.
        
        2. Among these candidate pockets, select the one that has:
           - The highest number of residues with contact changes
           - The highest number of interface residues
           - Good druggability scores (Score > 0.25 or Druggability Score > 0.05)
        
        3. For the selected pocket, include ALL residues from that pocket in your hotspots list, but order them by priority:
           - First: Residues that have both high contact changes AND are interface residues
           - Second: Residues that have either high contact changes OR are interface residues
           - Third: Remaining pocket residues
        
        The hotspots should be a comma-separated string of residue numbers with chain identifiers.
        ALL residues must come from the SAME pocket - do not mix residues from different pockets.
        The example residues are just for reference - use the actual residues from your analysis.
        
        IMPORTANT: Your selection MUST include:
        - ALL residues from the chosen pocket
        - At least one residue with high contact changes (total_change > 3.0)
        - At least one interface residue
        """
        
        print("  - Sending request to Anthropic API...")
        specifications_response = client.messages.create(
            model="claude-3-7-sonnet-20250219",
            max_tokens=1000,
            temperature=0.7,
            system="You are an expert in protein design and RFdiffusion.",
            messages=[{
                "role": "user",
                "content": specifications_prompt
            }]
        )
        print("  - Response received from Anthropic API")
        
        # Save API response to file
        api_output_file = os.path.join(OUTPUT_DIR, "anthropic_api_response.txt")
        with open(api_output_file, "w") as f:
            f.write("=== Anthropic API Response ===\n")
            f.write(f"Model: {specifications_response.model}\n")
            f.write(f"Usage: {specifications_response.usage}\n")
            f.write(f"Content: {specifications_response.content[0].text}\n")
        print(f"  - API response saved to {api_output_file}")
        
        content = specifications_response.content[0].text
        json_start = content.find("{")
        json_end = content.rfind("}") + 1
        
        if json_start == -1 or json_end == -1:
            raise ValueError("Could not find valid JSON in Anthropic response")
            
        return json.loads(content[json_start:json_end])
        
    except anthropic.APIError as e:
        print(f"⚠️ Anthropic API Error: {str(e)}")
        error_file = os.path.join(OUTPUT_DIR, "anthropic_api_error.txt")
        with open(error_file, "w") as f:
            f.write(f"=== Anthropic API Error ===\n")
            f.write(f"Error: {str(e)}\n")
        print(f"  - Error details saved to {error_file}")
        raise
    except json.JSONDecodeError as e:
        print(f"⚠️ Error parsing JSON from Anthropic response: {str(e)}")
        error_file = os.path.join(OUTPUT_DIR, "anthropic_api_error.txt")
        with open(error_file, "w") as f:
            f.write(f"=== JSON Decode Error ===\n")
            f.write(f"Error: {str(e)}\n")
        print(f"  - Error details saved to {error_file}")
        raise
    except Exception as e:
        print(f"⚠️ Unexpected error during specification generation: {str(e)}")
        error_file = os.path.join(OUTPUT_DIR, "anthropic_api_error.txt")
        with open(error_file, "w") as f:
            f.write(f"=== Unexpected Error ===\n")
            f.write(f"Error: {str(e)}\n")
        print(f"  - Error details saved to {error_file}")
        raise

def collect_specifications():
    """Collect specifications from the current analysis."""
    output = {}
    
    # Check for specifications.json in the current analysis directory
    spec_path = OUTPUT_DIR / "specifications.json"
    if spec_path.exists():
        try:
            with open(spec_path) as f:
                data = json.load(f)
            # Use S252W as the mutation name since that's what we're analyzing
            output["S252W"] = data.get("hotspots", "")
        except Exception as e:
            print(f"Warning: Could not read {spec_path}: {e}")

    # Create output file with datetime in analysis directory
    dt_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = OUTPUT_DIR / f"specifications_{dt_str}.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"✓ Combined specifications written to {out_path}")

def collect_contact_changes():
    """Collect contact changes from the current analysis."""
    output = {}
    
    # Check for most_changed_residues.json in the current analysis directory
    json_path = OUTPUT_DIR / "most_changed_residues.json"
    if json_path.exists():
        try:
            with open(json_path) as f:
                data = json.load(f)
            top5 = data.get("residues", [])[:5]
            output["S252W"] = top5
        except Exception as e:
            print(f"Warning: Could not read {json_path}: {e}")

    out_path = OUTPUT_DIR / "contact_changes.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"✓ Top 5 contact changes written to {out_path}")

def main():
    for pair in MUTATION_PAIRS:
        try:
            mut_path = pair["mut"]
            wt_path = pair["wt"]
            
            # Setup output directory
            output_dir = setup_output_dirs(mut_path)
            mutation_name = get_mutation_name(mut_path)
            protein_name = get_protein_name(wt_path)
            output_files = get_output_files(output_dir, protein_name, mutation_name)
            
            print(f"✓ Processing {mutation_name}...")
            
            # Load structures and create contact maps
            wt_struct = md.load(wt_path)
            mut_struct = md.load(mut_path)
            wt_map = ContactFrequency(wt_struct)
            mut_map = ContactFrequency(mut_struct)
            print("✓ Contact Maps Generated")
            
            # Create difference object and get contact differences
            diff = AtomMismatchedContactDifference(wt_map, mut_map)
            
            # Generate and save contact difference data
            top_positive = diff.residue_contacts.most_common()[:20]
            bottom_negative = list(reversed(diff.residue_contacts.most_common()))[:20]
            
            # Write contact differences to file
            with open(output_files["output_txt"], "w") as f:
                f.write("Top 20 positive differences:\n")
                for (pair, val) in top_positive:
                    f.write(f"{pair} => {val:.3f}\n")
                f.write("\nTop 20 negative differences:\n")
                for (pair, val) in bottom_negative:
                    f.write(f"{pair} => {val:.3f}\n")
            print("✓ Contact Differences File Generated")

            # Save difference matrix as CSV
            df = diff.residue_contacts.df
            df.to_csv(output_files["diff_csv"], index=True, header=True)
            print("✓ Contact Matrix CSV Generated")
            
            # Plot and save contact map
            fig, ax = diff.residue_contacts.plot()
            plt.title(f"Contact Map Difference: {protein_name} WT vs. {mutation_name}")
            plt.xlabel("Residue")
            plt.ylabel("Residue")
            plt.savefig(output_files["output_png"], dpi=150)
            plt.close()
            print("✓ Contact Map Plot Generated")

            # Analyze residue changes
            per_res_change = df.abs().sum(axis=1)
            changed_residues_series = per_res_change[per_res_change > 0].sort_values(ascending=False)

            # Convert indices to PDB labels
            changed_residues = []
            for idx in changed_residues_series.index:
                residue_obj = mut_struct.topology.residue(int(idx))
                
                chain_id = get_chain_id(residue_obj.chain)
                residue_str = f"{residue_obj.name} {chain_id}{residue_obj.resSeq}"
                
                changed_residues.append({
                    "rank": len(changed_residues) + 1,
                    "index": int(idx),
                    "pdb_label": residue_str,
                    "total_change": float(changed_residues_series[idx])
                })

            # Write contact analysis results to JSON
            output_data = {
                "analysis_type": "contact_map_hotspots",
                "description": (
                    "Residues in the mutant structure with any contact differences "
                    "compared to wild type. Based on non-zero values from the matrix row sums."
                ),
                "residues": changed_residues
            }

            with open(output_files["output_json"], "w") as f:
                json.dump(output_data, f, indent=2)
            print("✓ Contact Hotspots JSON Generated")

            # Run fpocket analysis on mutant structure
            mutant_pockets = analyze_pockets(mut_path)
            pocket_analysis = mutant_pockets
            print("✓ fpocket Analysis Complete")

            # Write pocket analysis to JSON
            with open(output_files["pocket_json"], "w") as f:
                json.dump(pocket_analysis, f, indent=2)
            print("✓ Pocket Analysis JSON Generated")

            # Get residue boundaries
            mutant_boundaries = get_residue_boundaries(mut_path)
            
            # Write to JSON file
            with open(output_files["boundaries_json"], "w") as f:
                json.dump(mutant_boundaries, f, indent=2)
            print("✓ Residue Boundaries JSON Generated")

            # Create combined analysis
            sequences = extract_sequences(wt_path, mut_path)
            
            # Save sequences separately
            sequences_file = output_files["sequences_json"]
            with open(sequences_file, "w") as f:
                json.dump(sequences, f, indent=2)
            print("✓ Sequences JSON Generated")
            
            # Combine the data
            combined_data = {
                "contact_analysis": output_data,
                "pocket_analysis": pocket_analysis,
                "residue_boundaries": mutant_boundaries,
                "sequences": sequences
            }
            
            # Write combined JSON
            with open(output_files["combined_json"], "w") as f:
                json.dump(combined_data, f, indent=2)
            print("✓ Combined Analysis JSON Generated")
            
            # Generate specifications
            hotspots = generate_specifications(combined_data)
            
            # Save specifications
            with open(output_files["spec_json"], "w") as f:
                json.dump(hotspots, f, indent=2)
            print("✓ Specifications JSON Generated")
            
            # Output summary
            base_name = Path(mut_path).stem
            fpocket_dir = output_dir / f"{base_name}_out"
            pocket_pdb = fpocket_dir / f"{base_name}_out.pdb"
            
            print(f"✓ Hotspots: {hotspots['hotspots']}")
            print(f"✓ PyMOL: {pocket_pdb}")
            print(f"✓ Output: {output_dir}")

        except Exception as e:
            mutation_name = get_mutation_name(mut_path)
            print(f"✗ Error: {mutation_name} - {str(e)}")
            
            # Save error details to file
            error_file = OUTPUT_DIR / "error_log.txt"
            with open(error_file, "w") as f:
                f.write(f"=== Error Processing {mutation_name} ===\n")
                f.write(f"Error: {str(e)}\n")
                f.write(f"Wildtype: {wt_path}\n")
                f.write(f"Mutant: {mut_path}\n")
            continue

if __name__ == "__main__":
    main()
    collect_specifications()
    collect_contact_changes()
