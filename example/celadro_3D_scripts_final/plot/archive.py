# This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import archive_base


class archive(archive_base.archive):
    """Simply reshape 2d fields after importing"""

    def __init__(self, path):
        super(archive, self).__init__(path)
        lx, ly, lz = self.Size
        self.parameters['walls'] = np.reshape(self.parameters['walls'],
                                              (lz,lx,ly))
        #self.parameters['walls'] = np.reshape(self.parameters['walls'],
        #                                     self.Size)
        self.__dict__.update(self.parameters)

    def read_frame(self, frame):
        frame = super(archive, self).read_frame(frame)

        # array sizes
        lx, ly, lz = self.Size
        px, py, pz = self.patch_size

        phi = []
        for i in range(len(frame.phi)):
            p = np.reshape(frame.phi[i], (pz, px, py))

            p = np.roll(p, frame.offset[i][0], axis=1)
            p = np.roll(p, frame.offset[i][1], axis=2)
            p = np.roll(p, frame.offset[i][2], axis=0)
            #print(p.shape)
            #print(lz-pz)
            #print(px)
            #print(py)
            p = np.concatenate((p, np.zeros((lz-pz, px,py))), axis=0)
            p = np.concatenate((p, np.zeros((lz,lx-px,py))), axis=1)
            p = np.concatenate((p, np.zeros((lz,lx,ly-py))), axis=2)      

            p = np.roll(p, frame.patch_min[i][0], axis=1)
            p = np.roll(p, frame.patch_min[i][1], axis=2)
            p = np.roll(p, frame.patch_min[i][2], axis=0)


            phi.append(p)
        frame.phi = phi

        if hasattr(frame, 'stress_xx'):
            frame.stress_xx.shape = (lz, lx, ly)
            frame.stress_xy.shape = (lz, lx, ly)
            frame.stress_yy.shape = (lz, lx, ly)
            frame.stress_zz.shape = (lz, lx, ly)
            frame.stress_xz.shape = (lz, lx, ly)

        return frame


def loadarchive(path):
    return archive(path)
