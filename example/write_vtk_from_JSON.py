import sys

sys.path.insert(0, "celadro_3D_scripts_final/plot/")
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
#import plot
import archive
#import animation
import gc

##################################################
# Init

if len(sys.argv) == 1:
    print("Please provide an input file.")
    exit(1)

# load archive from file
ar = archive.loadarchive(sys.argv[1])

oname = ""
if len(sys.argv) == 3:
    oname = "movie_"+sys.argv[2]
    print("Output name is", sys.argv[2])
    

def _get_field(phases, vals, size=1, mode='wrap'):
    """
    Compute the coarse grained field from a collection of phase-fields and
    associated values: ret[i] = sum_j phases[j]*values[i, j].

    Args:
        phases: List of phase-fields.
        vals: List of lists of size (None, len(phases)) of values to be
            associated with each phase-field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.

    Returns:
        A list of fields, each corresponding to the individual values.
    """
    ret = []
    for vlist in vals:
        assert len(vlist) == len(phases)
        field = np.zeros(phases[0].shape)
        for n in range(len(phases)):
            field += vlist[n]*phases[n]
        field = ndimage.filters.uniform_filter(field, size=size, mode=mode)
        ret.append(field)
    return ret


def get_velocity_field(phases, vel, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    v0 = [v[0] for v in vel]
    v1 = [v[1] for v in vel]
    v2 = [v[2] for v in vel]
    return _get_field(phases, [v0, v1, v2], size, mode)

    
def writeVTK(fname, frame):
    N = (frame.parameters['Size'][0]*frame.parameters['Size'][1]*frame.parameters['Size'][2])
    p = np.zeros(N)
    p = np.reshape(p,(frame.parameters['Size'][2],frame.parameters['Size'][0],frame.parameters['Size'][1]))
    for i in range(len(frame.phi)):
    	p += frame.phi[i]

    f = open(fname,"w+")
    f.write("# vtk DataFile Version 2.0\n");
    f.write("volume example\n");
    f.write("ASCII\n");
    f.write("DATASET STRUCTURED_POINTS\n");    
    f.write("DIMENSIONS %d %d %d\n" % (frame.parameters['Size'][1],frame.parameters['Size'][0],frame.parameters['Size'][2]));
    f.write("ASPECT_RATIO  %d %d %d\n" % (1, 1, 1));
    f.write("ORIGIN  %d %d %d\n" % (0, 0, 0));    
    f.write("POINT_DATA  %d\n" % N);
    f.write("SCALARS volume_scalars double 1\n");       
    f.write("LOOKUP_TABLE default \n");       

    for z in np.arange(0,frame.parameters['Size'][2],1):
    	for x in np.arange(0,frame.parameters['Size'][0],1):
    		for y in np.arange(0,frame.parameters['Size'][1],1):
    			f.write("%.6f\n" % p[z,x,y])

    del p
    gc.collect()
    f.close()
    
    
print('running write_vtk_from_JSON_05012021.py')
rng = np.arange(1,ar._nframes+1,1)
for fr in (rng):
    fname = 'frame_' + str(fr) + '.vtk'
    print(fname,flush=True)
    frame = ar.read_frame(fr)
    writeVTK(fname,frame)
    del frame
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
