import subprocess
import glob
import os

# you need steep.mdp in your working to run the grompp command
# I use a99SBdisp.ff to assign partial charges and hydrogens (in my working directory)
# you can probably run this with other force fields with minimal modifications
def get_charges_and_hydrogens(file, suppress=False):
    print(file)
    protein_ligand = file.split('_unrelaxed')[0]
    print(protein_ligand)

    # set file names
    mdp_file = "steep.mdp"
    conf_file = f"{protein_ligand}.gro"
    topol_file = f"{protein_ligand}.top"
    tpr = f"{protein_ligand}.tpr"
    output_pdb = f"{protein_ligand}.pdb"

    if suppress:
        stdout, stderr = subprocess.DEVNULL, subprocess.DEVNULL
    else:
        stdout, stderr = None, None

    try:
        # run pdb2gmx
        process = subprocess.Popen(
            ["gmx", "pdb2gmx", "-f", file, "-o",conf_file, "-p", topol_file],
            stdin=subprocess.PIPE,
            stdout=stdout,
            stderr=stderr,
            text=True
        )
        # Send the input and wait for the command to complete
        process.communicate(input="1\n1")

        # run editconf to center the protein in a large box to avoid periodic boundary issues
        process = subprocess.Popen(
            ["gmx", "editconf", "-f", conf_file, "-o", conf_file, "-c", "-d", "10", "-bt", "cubic"],
            stdin=subprocess.PIPE,
            stdout=stdout,
            stderr=stderr,
            text=True
        )
        process.communicate(input="\n")

        # run grompp to get a tpr file (contains charge info)
        process = subprocess.Popen(
            ["gmx", "grompp", "-f", mdp_file, "-c", conf_file, "-p", topol_file, "-o", tpr, "-maxwarn", "1"],
            stdin=subprocess.PIPE,
            stdout=stdout,
            stderr=stderr,
            text=True
        )
        process.communicate(input="\n")

        # convert back to pdb
        process = subprocess.Popen(
            ["gmx", "trjconv", "-f", conf_file, "-s", tpr, "-o", output_pdb],
            stdin=subprocess.PIPE,
            stdout=stdout,
            stderr=stderr,
            text=True
        )

        # Send the input '0' (followed by newline) to the process and wait for the command to complete
        process.communicate(input="0\n")

        # Delete specific files and those starting with #
        for pattern in ['#*', '*.top', 'posre*', 'mdout.mdp']:
            for file in glob.glob(pattern):
                os.remove(file)

        print(f"Outputted {output_pdb} and {tpr} (which contains charge info for hbond script)")
    except:
        print(f"Error processing {file}")

# Example usage
if __name__ == "__main__":
    suppress_output = False
    get_charges_and_hydrogens("HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb", suppress_output)
    # get_charges_and_hydrogens("pdbs/CO2PDZ3_COCRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb")
