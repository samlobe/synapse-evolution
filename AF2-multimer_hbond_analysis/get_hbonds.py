import pandas as pd
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def get_hbonds(file,debug=False):
    try:
        protein_ligand = file.split('_unrelaxed')[0]
        pdb = f'{protein_ligand}.pdb'
        tpr = f'{protein_ligand}.tpr'

        u = mda.Universe(tpr, pdb)

        hbonds = HBA(universe=u, between=['segid seg_0_Protein_chain_A', 'segid seg_1_Protein_chain_B'],
                        d_h_cutoff=1.2, d_a_cutoff=3.7, d_h_a_angle_cutoff=120)
        # add the backbone hydrogen
        hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
        hbonds.hydrogens_sel = f'{hbonds.hydrogens_sel} or (name H)'

        if debug:
            hbonds.run(verbose=True)
        else:
            hbonds.run()
        data = hbonds.results.hbonds

        # get resname, resid, chain, and atom name for the donor and hydrogen
        n_residues_A = len(u.select_atoms('segid seg_0_Protein_chain_A').residues)
        n_residues_B = len(u.select_atoms('segid seg_1_Protein_chain_B').residues)

        if debug:
            print(f"n_residues_A: {n_residues_A}")
            print(f"n_residues_B: {n_residues_B}")

        resnames_acceptor = [] ; resids_acceptor = [] ; chains_acceptor = []; names_acceptor = []
        resnames_hydrogen = [] ; resids_hydrogen = [] ; chains_hydrogen = []; names_hydrogen = []

        # parse the data
        for hbond in data:
            acceptor = int(hbond[3])
            hydrogen = int(hbond[2])
            names_acceptor.append(u.select_atoms(f'id {acceptor}').names[0])
            resnames_acceptor.append(u.select_atoms(f'id {acceptor}').resnames[0])
            resids_acceptor.append(u.select_atoms(f'id {acceptor}').resids[0])
            chains_acceptor.append(u.select_atoms(f'id {acceptor}').segids[0][-1])
            names_hydrogen.append(u.select_atoms(f'id {hydrogen}').names[0])
            resnames_hydrogen.append(u.select_atoms(f'id {hydrogen}').resnames[0])
            resids_hydrogen.append(u.select_atoms(f'id {hydrogen}').resids[0])
            chains_hydrogen.append(u.select_atoms(f'id {hydrogen}').segids[0][-1])

        # fix the resid numbering for chain B (gromacs issue)
        for i in range(len(resids_acceptor)):
            if chains_acceptor[i] == 'B':
                resids_acceptor[i] = resids_acceptor[i] - n_residues_A
        for i in range(len(resids_hydrogen)):
            if chains_hydrogen[i] == 'B':
                resids_hydrogen[i] = resids_hydrogen[i] - n_residues_A
        
        # get strings that descibe the hydrogen bonds for the last 6 ligand residues
        # rows: 0 site, -1 site, -2 site, -3 site, -4 site, -5 site
        # format should be something like "{one-letter-code}{PDZ-resid} {bb/sc}-{bb/sc}"
        # e.g. "H63 bb-sc" or  "I18 sc-sc"
                
        dict_3to1 = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
                    'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I',
                    'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
                    'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        
        if debug:
            print([f'{dict_3to1[resname]}{resid}{chain}' for resname, resid, chain in zip(resnames_acceptor, resids_acceptor, chains_acceptor)])
            print([f'{dict_3to1[resname]}{resid}{chain}' for resname, resid, chain in zip(resnames_hydrogen, resids_hydrogen, chains_hydrogen)])

        hbond_strings = []
        if debug:
            print(np.arange(n_residues_B+1)[-1:-7:-1])
        for resid in np.arange(n_residues_B+1)[-1:-7:-1]: # e.g. 14, 13, 12, 11, 10, 9
            res_hbonds_string = ''
            for name_a, resname_a, resid_a, chain_a, name_h, resname_h, resid_h, chain_h in zip(names_acceptor , resnames_acceptor, resids_acceptor, chains_acceptor, names_hydrogen, resnames_hydrogen, resids_hydrogen, chains_hydrogen):
                # if the acceptor is on chain B and hydrogen donor is chain A
                if chain_a == 'B' and resid_a == resid:
                    print(resid_a)
                    res_hbonds_string += f'{dict_3to1[resname_h]}{resid_h}'
                    if name_h == 'H': res_hbonds_string += 'bb-'
                    else: res_hbonds_string += 'sc-'

                    if name_a in ['O','OC1','OC2']: res_hbonds_string += 'bb' # considering C-terminus
                    else: res_hbonds_string += 'sc'
                    res_hbonds_string += ' '

                # if the hydrogen donor is on chain B and acceptor is chain A
                elif chain_h == 'B' and resid_h == resid:
                    res_hbonds_string += f'{dict_3to1[resname_a]}{resid_a}'
                    if name_a in ['O','OC1','OC2']: res_hbonds_string += 'bb-' # considering C-terminus
                    else: res_hbonds_string += 'sc-'

                    if name_h == 'H': res_hbonds_string += 'bb'
                    else: res_hbonds_string += 'sc'
                    res_hbonds_string += ' '

                else:
                    continue

            hbond_strings.append(res_hbonds_string)

        # create a dataframe
        hbond_df = pd.DataFrame(hbond_strings, columns=[f'{protein_ligand} hbonds'],index=[f'0 site', '-1 site', '-2 site', '-3 site', '-4 site', '-5 site'])
        hbond_df.index.name = 'ligand sites'
        return hbond_df
        # hbond_df.to_csv(f'hbonds_{protein_ligand}.csv')
    except:
        print(f"Error processing {file}")

if __name__ == "__main__":
    # hbond_df = get_hbonds("HSPDZ3_CRIPT_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb")
    # hbond_df = get_hbonds("pdbs/HSPDZ_CO_NEUROLIGIN_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb",debug=True)
    hbond_df = get_hbonds("pdbs/NODE11_HS_NEUROLIGIN_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb",debug=True)

    print(hbond_df)
