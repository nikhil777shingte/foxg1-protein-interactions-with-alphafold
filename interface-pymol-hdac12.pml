import json

# Define chains based on manual inspection or known data
hdac2_chain = 'B'
foxg1_chain = 'A' 

# Load the CIF file
cmd.load('/Users/nikhil/Documents/github/foxg1-hdac1-hdac2-jul2026/fold_foxg1_hdac2/fold_foxg1_hdac2_model_0.cif', 'model0')
#cmd.load('/Users/nikhil/Documents/github/foxg1-del-hdac2-2026/fold_181_275del_hfoxg1_hdac2/fold_181_275del_hfoxg1_hdac2_model_0.cif', 'model0')


# Select and color the HDAC2 protein
cmd.select('hdac2', f'chain {hdac2_chain} and model0')
cmd.show('cartoon', 'hdac2')
cmd.color('blue', 'hdac2')

# Select and color the FOXG1 protein
cmd.select('foxg1', f'chain {foxg1_chain} and model0')
cmd.show('cartoon', 'foxg1')
cmd.color('red', 'foxg1')

# Select and show interacting residues within 4 Å
interface = interfaceResidues('model0', 'c. A', 'c. B', cutoff=4.0, selName="interface")
cmd.show('sticks', 'interface')

# Select and color the interface residues for HDAC2
cmd.select('hdac2_interface', f'interface and chain {hdac2_chain}')
cmd.color('blue', 'hdac2_interface')

# Select and color the interface residues for FOXG1
cmd.select('foxg1_interface', f'interface and chain {foxg1_chain}')
cmd.color('red', 'foxg1_interface')

cmd.zoom('interface')

# Get the list of interface residues for HDAC2
hdac2_residues = set()
cmd.iterate('hdac2_interface and not elem H', 'hdac2_residues.add((resn, resi))')

# Get the list of interface residues for FOXG1
foxg1_residues = set()
cmd.iterate('foxg1_interface and not elem H', 'foxg1_residues.add((resn, resi))')

# Convert sets to sorted lists for readability
hdac2_residues = sorted(hdac2_residues, key=lambda x: int(x[1]))
foxg1_residues = sorted(foxg1_residues, key=lambda x: int(x[1]))

# Print the interface residues for each protein
print("HDAC2 interface residues:")
print([f"{resn}{resi}" for resn, resi in hdac2_residues])

print("FOXG1 interface residues:")
print([f"{resn}{resi}" for resn, resi in foxg1_residues])
