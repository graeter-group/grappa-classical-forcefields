#%%
from kimmdy.topology.topology import Topology
from kimmdy.parsing import read_top
from pathlib import Path
import math
import numpy as np

#%%
print("Converting GROMACS ff directory to OpenMM xml file")
ff_type = 'amber'  # charmm or amber
stem = "grappa_1-3-amber99_unique"
#stem = "grappa_1-3-amber99_ff19SB_trimersL"
#stem = 'amber99sb-star-ildnp_mod'
ffin = Path(f"{stem}.ff")
ffout = Path(f"{stem}.xml")
ATNUM_ELE = {'1':'H','6':'C','7':'N','8':'O','9':'F','35':'Br','20':'Ca','53':'I','17':'Cl','11':'Na','0':'Du','12':'Mg','15':'P','16':'S','29':'Cu','26':'Fe','19':'K','37':'Rb','55':'Cs','3':'Li','30':'Zn'}
print(f"GROMACS ff directory: {ffin} of type {ff_type}, writing to {ffout}.")

#%%
# create top file with includes to desired force field
top_template = Path('top_templates/generic_template.top')
top_ff = Path(f"top_templates/{stem}.top")
print(f"Writing template file for KIMMDY parsing at {top_ff}")

with open(top_template, 'r') as file:
    filedata = file.read()

filedata = filedata.replace('XX', ffin.name)

with open(top_ff, 'w') as file:
    file.write(filedata)

#%%
print("Reading GROMACS topology file using KIMMDY")
ff = Topology(read_top(top_ff,ffdir=ffin),residuetypes_path=ffin / 'aminoacids.rtp',radicals=' ') #take any topology, doesn't matter, change residuetypes_path if necessary

#%%
print('Extracting information for force field construction')
## reduce residues to protein part
protein_residues = ["ACE","NME","ALA", "ARG", "ASN", "ASP", "ASPP","ASH","CYS", "CYN","CYM","GLN", "GLU","GLH","GLUP", "GLY","HYP",
                     "HSE","HSP","HSD", "HID","HIP","HIE","ILE","LEU", "LYS","LSN","LYN","MET","PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]#,"DOP"
# ,"CYX" # is disulfide, CYM is deprotonated
protein_residues.extend(['C' + res for res in protein_residues])
protein_residues.extend(['N' + res for res in protein_residues])
# "NALA", "NARG", "NASN", "NASP","NASPP", "NCYS", "NGLN", "NGLU", "NGLUP","NGLY", "NHSE", "NHSP", "NHSD", "NILE", "NLEU", "NLYS", "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL","CALA", "CARG", "CASN", "CASP","CASPP", "CCYS", "CGLN", "CGLU", "CGLUP","CGLY", "CHSE", "CHSP", "CHSD", "CILE","CLEU", "CLYS", "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP",  "CTYR", "CVAL"
atomtypes_protein = set()

rmv_res = []
for residue,entry in ff.ff.residuetypes.items():
    if residue not in protein_residues:
        rmv_res.append(residue)
print(f"Removing residues {rmv_res}; these won't be part of the xml file!")
for res in rmv_res:
    del ff.ff.residuetypes[res]


for residue,entry in ff.ff.residuetypes.items():
    if residue in protein_residues:
        atomtypes_protein |= set([atom.type for atom in entry.atoms.values()])
atomtypes_protein |= set(['1','2','3','4','5','6','7','8','9','0','X'])

#%%
interactiontype_dicts = [ff.ff.bondtypes,ff.ff.angletypes,ff.ff.proper_dihedraltypes,ff.ff.improper_dihedraltypes]
interactions = ['bond','angle','proper','improper']
for i,interactiontype_dict in enumerate(interactiontype_dicts):
    print(f"# {interactions[i]}types before removal: {len(interactiontype_dict)}")
    rmv_interactiontypes = []
    for interactionstype in interactiontype_dict.keys():
        if not all(x in atomtypes_protein for x in interactionstype):
            rmv_interactiontypes.append(interactionstype)
    for key in rmv_interactiontypes:
            del interactiontype_dict[key]
    print(f"# {interactions[i]}types after removal: {len(interactiontype_dict)}")

# %%
Atom = []
Residue = {}
i = 0
for rt_name, rt in ff.ff.residuetypes.items():
    Residue[rt_name] = {'Atom':[],'Bond':[],'ExternalBond':[]}
    residue_name_to_idx = {}
    for j, (atom, atomspec) in enumerate(rt.atoms.items()):
        at = atomspec.type
        type_info = ff.ff.atomtypes[at]
        Atom.append({'name':str(i),'class':at,'element':ATNUM_ELE[type_info.at_num],'mass':type_info.mass,'charge':atomspec.charge,'sigma':type_info.sigma,'epsilon':type_info.epsilon})
        Residue[rt_name]['Atom'].append({'name':atomspec.name,'type':str(i)})
        residue_name_to_idx[atomspec.name] = str(j)
        if atomspec.name == 'C':
            if not (len(rt_name) == 4 and rt_name.startswith('C')):
                Residue[rt_name]['ExternalBond'].append({'from':str(j)})
        if atomspec.name == 'N':
            if not (len(rt_name) == 4 and rt_name.startswith('N')):
                Residue[rt_name]['ExternalBond'].append({'from':str(j)})
        i += 1
    for bond, bondspec in rt.bonds.items():
        # this means the bond is between different residues (before/after)
        if any([x[0] in ['+','-'] for x in bond]):
            fromatom = bond[0] if bond[0][0] not in ['+','-'] else bond[1]
            if not any ([x['from'] == residue_name_to_idx[fromatom] for x in Residue[rt_name]['ExternalBond']]):
                Residue[rt_name]['ExternalBond'].append({'from':residue_name_to_idx[fromatom]})
        else:
            Residue[rt_name]['Bond'].append({'from':residue_name_to_idx[bond[0]],'to':residue_name_to_idx[bond[1]]})

#%%
## charmm specific ##
if ff_type == 'charmm':
    ## do cmaps separately until kimmdy can deal with them
    with open(ffin/"cmap.itp",'r') as f:
        foo = f.readlines()
    cmap = {}
    for line in foo:
        linestrip = line.strip()
        linestrip = linestrip.strip('\\')
        if linestrip.startswith(';') or linestrip.startswith('[') or len(linestrip) == 0:
            continue
        linesplit = linestrip.split()
        if len(linesplit) == 8:
            key = tuple(linesplit[:5])
            cmap[key] = []
        else:
            cmap[key].extend([str(x) for x in linesplit])

    cmap_shaped = {}
    for atoms, vals in cmap.items():
        cmap_shaped[atoms] = np.asarray(vals).reshape(24,24)

# %%
## write to xml file ##
print("Writing xml output.")
with open(ffout,'w') as f:
    f.write("<ForceField>\n")

    f.write(" <AtomTypes>\n")
    for atom in Atom:
        f.write(f"  <Type name=\"{atom['name']}\" class=\"{atom['class']}\" element=\"{atom['element']}\" mass=\"{atom['mass']}\"/>\n")
    f.write(" </AtomTypes>\n")

    f.write(" <Residues>\n")
    for resname, res in Residue.items():
        f.write(f"  <Residue name=\"{resname}\">\n")
        for atom in res['Atom']:
            f.write(f"   <Atom name=\"{atom['name']}\" type=\"{atom['type']}\"/>\n")
        for bond in res['Bond']:
            f.write(f"   <Bond from=\"{bond['from']}\" to=\"{bond['to']}\"/>\n")
        for ext_bond in res['ExternalBond']:
            f.write(f"   <ExternalBond from=\"{ext_bond['from']}\"/>\n")
        f.write(f"  </Residue>\n")
    f.write(" </Residues>\n")

    f.write(f" <HarmonicBondForce>\n")
    for interactionstype in ff.ff.bondtypes.values():
        f.write(f"  <Bond class1=\"{interactionstype.i}\" class2=\"{interactionstype.j}\" length=\"{interactionstype.c0}\" k=\"{interactionstype.c1}\"/>\n")
    f.write(f" </HarmonicBondForce>\n")

    f.write(f" <HarmonicAngleForce>\n")
    for angletype in ff.ff.angletypes.values():
        f.write(f"  <Angle class1=\"{angletype.i}\" class2=\"{angletype.j}\" class3=\"{angletype.k}\" angle=\"{float(angletype.c0)* (math.pi/180)}\" k=\"{angletype.c1}\"/>\n")
    f.write(f" </HarmonicAngleForce>\n")

    f.write(f" <AmoebaUreyBradleyForce>\n")
    for angletype in ff.ff.angletypes.values():
        if angletype.funct == '5':
            f.write(f"  <UreyBradley class1=\"{angletype.i}\" class2=\"{angletype.j}\" class3=\"{angletype.k}\" angle=\"{float(angletype.c2)* (math.pi/180)}\" k=\"{angletype.c3}\"/>\n")
    f.write(f" </AmoebaUreyBradleyForce>\n")

    if ff_type == 'amber':
        f.write(f" <PeriodicTorsionForce ordering = \"amber\">\n")
    else:
        f.write(f" <PeriodicTorsionForce>\n")
    dihedrals_by_atoms = {}
    for k,proper_dihedraltype in ff.ff.proper_dihedraltypes.items():
        k_atoms = k[:4]
        k_peridicity = k[4]
        if k_atoms not in dihedrals_by_atoms:
            dihedrals_by_atoms[k_atoms] = {}
        dihedrals_by_atoms[k_atoms][k_peridicity] = proper_dihedraltype

    for k,v in dihedrals_by_atoms.items():
        #i = '' if proper_dihedraltype.i == "X" else proper_dihedraltype.i
        k = ['' if x == "X" else x for x in k]
        f.write(f"  <Proper class1=\"{k[0]}\" class2=\"{k[1]}\" class3=\"{k[2]}\" class4=\"{k[3]}\"")
        for n, key in enumerate(v.keys()):
            f.write(f" periodicity{n+1}=\"{v[key].periodicity}\" phase{n+1}=\"{float(v[key].c0) * (math.pi/180)}\" k{n+1}=\"{v[key].c1}\"")
        f.write("/>\n")
    
    for improper_dihedraltype in ff.ff.improper_dihedraltypes.values():
        if improper_dihedraltype.funct == '4':
            i = '' if improper_dihedraltype.i == "X" else improper_dihedraltype.i
            j = '' if improper_dihedraltype.j == "X" else improper_dihedraltype.j
            k = '' if improper_dihedraltype.k == "X" else improper_dihedraltype.k
            l = '' if improper_dihedraltype.l == "X" else improper_dihedraltype.l
            f.write(f"  <Improper class1=\"{k}\" class2=\"{j}\" class3=\"{i}\" class4=\"{l}\" periodicity1=\"{improper_dihedraltype.periodicity}\" phase1=\"{float(improper_dihedraltype.c0) * (math.pi/180)}\" k1=\"{improper_dihedraltype.c1}\"/>\n")
    f.write(f" </PeriodicTorsionForce>\n")

    if ff_type == 'charmm':
        f.write(f" <CustomTorsionForce energy=\"k*(theta-theta0)^2\">\n  <PerTorsionParameter name=\"k\"/>\n  <PerTorsionParameter name=\"theta0\"/>\n")
        for improper_dihedraltype in ff.ff.improper_dihedraltypes.values():
            if improper_dihedraltype.funct == '2':
                i = '' if improper_dihedraltype.i == "X" else improper_dihedraltype.i
                j = '' if improper_dihedraltype.j == "X" else improper_dihedraltype.j
                k = '' if improper_dihedraltype.k == "X" else improper_dihedraltype.k
                l = '' if improper_dihedraltype.l == "X" else improper_dihedraltype.l
                f.write(f"  <Improper class1=\"{i}\" class2=\"{j}\" class3=\"{k}\" class4=\"{l}\" k=\"{improper_dihedraltype.c1}\" theta0=\"{improper_dihedraltype.c0}\"/>\n")
        f.write(f" </CustomTorsionForce>\n")
    
        f.write(f" <CMAPTorsionForce>\n")
        for vals in cmap_shaped.values():
            f.write("   <Map>")
            for row in vals:
                f.write(" " +" ".join(row) + "\n")
            f.write("</Map>\n")
        for i,keys in enumerate(cmap_shaped.keys()):
            f.write(f"    <Torsion map=\"{i}\" class1=\"{keys[0]}\" class2=\"{keys}\" class3=\"{keys[2]}\" class4=\"{keys[3]}\" class5=\"{keys[4]}\"/>\n")
        f.write(f" </CMAPTorsionForce>\n")

    f.write(f" <NonbondedForce coulomb14scale=\"0.833333\" lj14scale=\"0.5\">\n")
    for atom in Atom:
        f.write(f"  <Atom type=\"{atom['name']}\" charge=\"{atom['charge']}\" sigma=\"{atom['sigma']}\" epsilon=\"{atom['epsilon']}\"/>\n")
    f.write(f" </NonbondedForce>\n")
    f.write(f"</ForceField>")
print('Done!')
# %%
