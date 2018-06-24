import MDAnalysis as mda
import numpy as np

u = mda.Universe.empty(9, 3, atom_resindex=[0, 0, 0, 0, 1, 1, 2, 2, 2], trajectory=True)

# u.write('./test/new.data')

# with mda.Writer('./test/new.data', n_atoms=u.atoms.n_atoms) as W:
#     W.write(ts)
#     # for ts in u.trajectory:
#     #     W.write(ts)


u.add_TopologyAttr('names')
u.add_TopologyAttr('masses')
u.add_TopologyAttr('charges')
u.add_TopologyAttr('types')
Type_array = np.array(['1']*9)
u.atoms.types = Type_array

mysel = u.select_atoms("all")
mysel.write('./test/new.data')



N_frames = 1000
coordinates = np.ones((N_frames, u.atoms.n_atoms, 3))

u.load_new(coordinates, order="fac")

print(u.trajectory.n_frames)

# Now you have an empty trajectory attached that you can fill with your data. Or have coordinates prepared with your data.

# Write to a new format:

with mda.Writer('./test/new.dcd', n_atoms=u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(ts)
