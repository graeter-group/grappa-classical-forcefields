#%%
import shutil
import subprocess as sp
from pathlib import Path
from kimmdy.topology.topology import Topology
from kimmdy.topology.atomic import ResidueImproperSpec
from kimmdy.parsing import read_top
from copy import deepcopy

#%% constants
at_num = {'C':6,'O':8,'N':7,'H':1,'S':16}
mass = {'C':12.01,'O':16.00,'N':14.01,'H':1.008,'S':16.00}
AA = {
    'A': 'ALA',
    'C': 'CYS',
    'D': 'ASP',
    'E': 'GLU',
    'F': 'PHE',
    'G': 'GLY',
    'H': 'HIE',
    'I': 'ILE',
    'K': 'LYS',
    'L': 'LEU',
    'M': 'MET',
    'N': 'ASN',
    'P': 'PRO',
    'Q': 'GLN',
    'R': 'ARG',
    'S': 'SER',
    'T': 'THR',
    'V': 'VAL',
    'W': 'TRP',
    'Y': 'TYR',
    'J': 'DOP',
    'O': 'HYP',
    '1': 'HID',
    '2': 'HIP',
    '3': 'GLH',
    '4': 'ASH',
    '5': 'LYN',
    '6': 'CYM',
    '7': 'LEU',
    'B': 'ACE',
    'Z': 'NME',
    '': '',
}
AA3_remap = {
    'ALA': 'ALA',
    'CYS': 'CYS',
    'ASP': 'ASP',
    'GLU': 'GLU',
    'PHE': 'PHE',
    'GLY': 'GLY',
    'HISE': 'HIE',
    'ILE': 'ILE',
    'LYS': 'LYS',
    'LEU': 'LEU',
    'MET': 'MET',
    'ASN': 'ASN',
    'PRO': 'PRO',
    'GLN': 'GLN',
    'ARG': 'ARG',
    'SER': 'SER',
    'THR': 'THR',
    'VAL': 'VAL',
    'TRP': 'TRP',
    'TYR': 'TYR',
    'DOP': 'DOP',
    'HYP': 'HYP',
    'HISD': 'HID',
    'HISH': 'HIP',
    'GLUH': 'GLH',
    'ASPH': 'ASH',
    'LYSN': 'LYN',
    'CYM': 'CYM',
    'ACE': 'ACE',
    'NME': 'NME'
}
skip_impropers = False

#%% functions
def order_residues(res:str):
    if len(res) == 3:
        return (0, res)
    elif len(res) == 4 and res.startswith('N'):
        return (1, res[1:])  # Remove the first character 'N' for sorting
    elif len(res) == 4 and res.startswith('C'):
        return (2, res[1:])  # Remove the first character 'C' for sorting
    else:
        return (3, res)

def concat_files(src1:str,src2:str,dest:str,cwd:Path):
    shutil.copy2(cwd / src1,cwd / dest)
    sp.run(["cat", src2], stdout=open(cwd / dest, 'ab'),cwd=cwd)

def read_yml(filepath: Path):
    with open(filepath,'r') as f:
        foo = f.readlines()
    yml_dict = {}
    for line in foo:
        linesplit = [x .strip() for x in line.split(sep=":")]
        yml_dict[linesplit[0]] = linesplit[1]
    return yml_dict

#%% hyperparameters
top_dir = Path("topologies_grappa-1.3.0_amber99")
ff_dir = Path(f"grappa_1-3-amber99_unique.ff")
neighbor_AA = 'L'

datasets = ["sidechain","nter","cter","caps","ss-interaction","ns-interaction","cs-interaction"]   # always go through single dataset first, then interaction dataset
#shutil.copytree('grappa_blank_unique.ff',ff_dir,dirs_exist_ok=False)
config_file = 'dataset.yml'

#%%
atomtypes = {}
bondedtypes = {'bondtypes':{},'angletypes':{},'dihedraltypes_proper':{},'dihedraltypes_improper':{}}
residuetypes = {}

#%%
for dataset in datasets:
    config = read_yml(top_dir / dataset / config_file)
    top_paths = list((top_dir / dataset).glob("*top"))

    for top_str in top_paths:
        seq = top_str.stem
        ele_counter = {'C':1,'H':1,'N':1,'O':1,'S':1}
        if config['resid'] == 'cap':
            if seq == 'B':
                resid = ['1']
            elif seq == 'Z':
                resid = ['4']
            else:
                raise ValueError
        else:
            resid = list(config['resid'].split(sep=','))
        if config['res'] == 'seq':
            res = seq
        elif config['res'] == 'None':
            res = ''
        else:
            raise ValueError(f"Unexpected res value in config: {config['res']}")
        if config['resname_prefix'] != 'None':
            prefix = config['resname_prefix']
        else:
            prefix = ''
        resname = prefix + AA[res]

        if config['type'] not in ['single','interaction']:
            raise ValueError

        print(top_str, res, resname, resid)

        #read top
        top = Topology(read_top(top_str),residuetypes_path=Path("amber99sb-star-ildnp_mod.ff/aminoacids.rtp"))

        # change atom type for single amino acid parametrization
        if config['type'] == 'single':
            for type_res in top.atoms.values():
                if type_res.resnr == str(resid[0]):
                    # get sigma, epsilon, at_num, mass from old atomtype
                    type_info =top.ff.atomtypes[type_res.type]
                    info_raw = deepcopy({'charge':type_info.charge,'epsilon':type_info.epsilon,'sigma':type_info.sigma})
                    # construct new atomtype
                    ele = type_res.type[0]
                    atomtype = f"{ele}{res.lower()}{prefix.lower()}{ele_counter[ele]}"
                    ele_counter[ele] +=1
                    type_res.type = atomtype
                    #print(atomtype)
                    # fill dict for new ffnonbonded.itp
                    atomtypes[atomtype] = {'at_num':type_info.at_num,'mass':type_info.mass,'charge':'0.0000','ptype':'A','sigma':type_info.sigma,'epsilon':type_info.epsilon}
                    # fill dict for new aminoacids.rtp
                    top.ff.residuetypes[resname].atoms[type_res.atom].type = atomtype

        # bond part for ffbonded.itp
        for bond in top.bonds.values():
            atoms_info = {'type':[],'name':[],'resid':[],'resname':[]}
            for nr in [bond.ai,bond.aj]:
                atoms_info['type'].append(top.atoms[nr].type)
                atoms_info['name'].append(top.atoms[nr].atom)
                atoms_info['resid'].append(top.atoms[nr].resnr)
                atoms_info['resname'].append(top.atoms[nr].residue)
            if all([x in resid for x in atoms_info['resid']]):
                if config['type'] == 'single':
                    atoms_types = atoms_info['type']
                elif config['type'] == 'interaction':
                    if any([x == resid[0] for x in atoms_info['resid']]) and any([x == resid[1] for x in atoms_info['resid']]):
                        atoms_types = [residuetypes[AA3_remap[atoms_info['resname'][i]]].atoms[x].type for i,x in enumerate(atoms_info['name'])]
                    else:
                        # take parameters for interaction topology only if both residues are involved
                        continue
                bondedtypes['bondtypes'][tuple(atoms_types)] = {'funct':bond.funct,'c0':bond.c0,'c1':bond.c1}

        # angles
        for angle in top.angles.values():
            atoms_info = {'type':[],'name':[],'resid':[],'resname':[]}
            for nr in [angle.ai,angle.aj,angle.ak]:
                atoms_info['type'].append(top.atoms[nr].type)
                atoms_info['name'].append(top.atoms[nr].atom)
                atoms_info['resid'].append(top.atoms[nr].resnr)
                atoms_info['resname'].append(top.atoms[nr].residue)
            if all([x in resid for x in atoms_info['resid']]):
                if config['type'] == 'single':
                    atoms_types = atoms_info['type']
                elif config['type'] == 'interaction':
                    if any([x == resid[0] for x in atoms_info['resid']]) and any([x == resid[1] for x in atoms_info['resid']]):
                        atoms_types = [residuetypes[AA3_remap[atoms_info['resname'][i]]].atoms[x].type for i,x in enumerate(atoms_info['name'])]
                    else:
                        # take parameters for interaction topology only if both residues are involved
                        continue
                bondedtypes['angletypes'][tuple(atoms_types)] = {'funct':angle.funct,'c0':angle.c0,'c1':angle.c1}

        # proper dihedrals
        for multiple_dihedrals in top.proper_dihedrals.values():
            atoms_info = {'type':[],'name':[],'resid':[],'resname':[]}
            for nr in [multiple_dihedrals.ai,multiple_dihedrals.aj,multiple_dihedrals.ak,multiple_dihedrals.al]:
                atoms_info['type'].append(top.atoms[nr].type)
                atoms_info['name'].append(top.atoms[nr].atom)
                atoms_info['resid'].append(top.atoms[nr].resnr)
                atoms_info['resname'].append(top.atoms[nr].residue)
            if all([x in resid for x in atoms_info['resid']]):
                if config['type'] == 'single':
                    atoms_types = atoms_info['type']
                elif config['type'] == 'interaction':
                    if any([x == resid[0] for x in atoms_info['resid']]) and any([x == resid[1] for x in atoms_info['resid']]):
                        atoms_types = [residuetypes[AA3_remap[atoms_info['resname'][i]]].atoms[x].type for i,x in enumerate(atoms_info['name'])]
                    else:
                        # take parameters for interaction topology only if both residues are involved
                        continue
                for dihedral in multiple_dihedrals.dihedrals.values():
                    bondedtypes['dihedraltypes_proper'][tuple(atoms_types)+(dihedral.periodicity,)] = {'funct':dihedral.funct,'c0':dihedral.c0,'c1':dihedral.c1,'periodicity':dihedral.periodicity}
        # improper dihedrals
        for multiple_dihedrals in top.improper_dihedrals.values():
            atoms_info = {'type':[],'name':[],'resid':[],'resname':[]}
            for nr in [multiple_dihedrals.ai,multiple_dihedrals.aj,multiple_dihedrals.ak,multiple_dihedrals.al]:
                atoms_info['type'].append(top.atoms[nr].type)
                atoms_info['name'].append(top.atoms[nr].atom)
                atoms_info['resid'].append(top.atoms[nr].resnr)
                atoms_info['resname'].append(top.atoms[nr].residue)
            if all([x in resid for x in atoms_info['resid']]):
                if config['type'] == 'single':
                    atoms_types = atoms_info['type']
                elif config['type'] == 'interaction':
                    if any([x == resid[0] for x in atoms_info['resid']]) and any([x == resid[1] for x in atoms_info['resid']]):
                        atoms_types = [residuetypes[AA3_remap[atoms_info['resname'][i]]].atoms[x].type for i,x in enumerate(atoms_info['name'])]
                    else:
                        # take parameters for interaction topology only if both residues are involved
                        continue
                for dihedral in multiple_dihedrals.dihedrals.values():
                    bondedtypes['dihedraltypes_improper'][tuple(atoms_types)+(dihedral.periodicity,)] = {'funct':dihedral.funct,'c0':dihedral.c0,'c1':dihedral.c1,'periodicity':dihedral.periodicity}
                    # also add mirrored improper for cases where Gromacs switches impropers around (should gromacs be able to mirror by itself?)
                    #bondedtypes['dihedraltypes_improper'][tuple(atoms_types[::-1])+(dihedral.periodicity,)] = {'funct':dihedral.funct,'c0':dihedral.c0,'c1':dihedral.c1,'periodicity':dihedral.periodicity}
                    # and permutations (only switches sign)
                    bondedtypes['dihedraltypes_improper'][tuple([atoms_types[k] for k in [0,2,1,3]])+(dihedral.periodicity,)] = {'funct':dihedral.funct,'c0':dihedral.c0,'c1':dihedral.c1,'periodicity':dihedral.periodicity}
                    #bondedtypes['dihedraltypes_improper'][tuple([atoms_types[k] for k in [3,2,1,0]])+(dihedral.periodicity,)] = {'funct':dihedral.funct,'c0':dihedral.c0,'c1':dihedral.c1,'periodicity':dihedral.periodicity}

        if config['type'] == 'single':
            # add residuetype with new atomtypes to dict
            del_resimproper = []
            for k, improper in top.ff.residuetypes[resname].improper_dihedrals.items():
                improper.c0 = None
                improper.c1 = None
                improper.c2 = None
                if k ==  ('N','CA','C','+N'):
                    del_resimproper.append(k)   # linear improper
            for resimproper in del_resimproper:
                #print('Removed residue improper:')
                improper_popped = top.ff.residuetypes[resname].improper_dihedrals.pop(resimproper)
                # print(improper_popped)

            ## could also take actual impropers of the topology and write them as residueimproperspec instead
            add_resimproper = {}
            for k, improper in top.ff.residuetypes[resname].improper_dihedrals.items():
                #print(k)
                for permutation in [[1,3,2,0],[3,0,2,1]]:
                    k_new = tuple([k[p] for p in permutation])
                    add_resimproper[k_new] = ResidueImproperSpec.from_top_line(list(k_new))
                    #print(k_new)
            #print(add_resimproper)
            top.ff.residuetypes[resname].improper_dihedrals.update(add_resimproper)
            #print(list(top.ff.residuetypes[resname].improper_dihedrals.keys()))

            residuetypes[resname] = top.ff.residuetypes[resname]

# %%
with open(ff_dir / "ffnonbonded.itp",'w') as f:
    f.write("[ atomtypes ]\n")
    f.write("; name      at.num  mass     charge ptype  sigma      epsilon\n")
    for k,v in sorted(atomtypes.items()):
        f.write("{0:<5}        {1:<3}    {2:<6}   {3:<6}  {4}   {5:<10}   {6:<10}\n".format(k,*list(v.values())))
    f.write('#include "ffnonbonded_basic.itp"')

with open(ff_dir / "atomtypes_grappa.atp",'w') as f:
    for k,v in sorted(atomtypes.items()):
        f.write("{0:<5}             {1:<8}  \n".format(k,v['mass']))
#concat_files('atomtypes_basic.atp',"atomtypes_grappa.atp","atomtypes.atp",ff_dir)  #gromacs finds both files, no need to cat

with open(ff_dir / "ffbonded_bondtypes.itp",'w') as f:
    for k,v in sorted(bondedtypes['bondtypes'].items(),key=None):
        f.write(f"  {k[0]:>4} {k[1]:>4}        {v['funct']}  {float(v['c0']):>8.6f} {float(v['c1']):>10.2f}\n")

with open(ff_dir / "ffbonded_angletypes.itp",'w') as f:
    for k,v in sorted(bondedtypes['angletypes'].items(),key=None):
        f.write(f"  {k[0]:>4} {k[1]:>4} {k[2]:>4}       {v['funct']}  {float(v['c0']):>7.3f} {float(v['c1']):>9.3f}\n")

with open(ff_dir / "ffbonded_dihedraltypes_proper.itp",'w') as f:
    for k,v in sorted(bondedtypes['dihedraltypes_proper'].items(),key=None):
        f.write(f"  {k[0]:>4} {k[1]:>4} {k[2]:>4} {k[3]:>4}   {v['funct']}  {float(v['c0']):>7.3f} {float(v['c1']):>9.5f}   {v['periodicity']}\n")

with open(ff_dir / "ffbonded_dihedraltypes_improper.itp",'w') as f:
    for k,v in sorted(bondedtypes['dihedraltypes_improper'].items(),key=None):
        # tabulated force fields in gromacs only work with one improper per set of atoms, use periodicity = 2 because it is physically meaningful
        if skip_impropers and v['periodicity'] != '2':
            continue
        elif v['periodicity'] != '2':
            continue
        f.write(f"  {k[0]:>4} {k[1]:>4} {k[2]:>4} {k[3]:>4}   {v['funct']}  {float(v['c0']):>7.3f} {float(v['c1']):>9.5f}   {v['periodicity']}\n")

with open(ff_dir / "aminoacids.rtp",'w') as f:
    f.write("[ bondedtypes ]\n")
    f.write("; Col 1: Type of bond\n")
    f.write("; Col 2: Type of angles\n")
    f.write("; Col 3: Type of proper dihedrals\n")
    f.write("; Col 4: Type of improper dihedrals\n")
    f.write("; Col 5: Generate all dihedrals if 1, only heavy atoms of 0.\n")
    f.write("; Col 6: Number of excluded neighbors for nonbonded interactions\n")
    f.write("; Col 7: Generate 1,4 interactions between pairs of hydrogens if 1\n")
    f.write("; Col 8: Remove impropers over the same bond as a proper if it is 1\n")
    f.write("; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih\n")
    f.write("     1       1          9          4        1         3      1     0\n")
    f.write("\n")
    f.write("; now: water, ions, urea, terminal caps, AA's and terminal AA's\n")
    f.write("\n")
    for k,v in sorted(residuetypes.items(),key=lambda x: order_residues(x[0])):
        print(k)
        f.write(f"[ {k} ]\n")
        #atoms
        f.write(f" [ atoms ]\n")
        for type_res in v.atoms.values():
            f.write(f" {type_res.name:>4} {type_res.type:>5}          {type_res.charge:>9}  {type_res.cgrp:>2}\n")
        #bonds
        f.write(f" [ bonds ]\n")
        for bond in v.bonds.values():
            f.write(f" {bond.atom1:>4}   {bond.atom2:>4}\n")
        #impropers
        f.write(f" [ impropers ]\n")
        for improper in v.improper_dihedrals.values():
            f.write(f" {improper.atom1:>4}  {improper.atom2:>4}  {improper.atom3:>4}  {improper.atom4:>4}")
            for val in [improper.c0,improper.c1,improper.c2]:
                if val is not None:
                    f.write(f" {val:>7}")
            f.write("\n")
        f.write("\n")
# concat_files('aminoacids_grappa.rtp',"aminoacids_basic.rtp","aminoacids.rtp",ff_dir) #gromacs finds both files, no need to cat
# %%
