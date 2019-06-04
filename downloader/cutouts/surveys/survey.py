import os
import tempfile
import shutil

import montage_wrapper

# abstract class for a survey
from abc import ABC, abstractmethod
class Survey(ABC):

    def __init__(self):
        super().__init__()

    # grab a cutout of size <size> centered on <position>
    @abstractmethod
    def get_cutout(self, position, size):
        pass

    # make the directory structure if it doesn't exist
    def __make_dir(self, dirname):
        try:
            os.makedirs(dirname)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    def mosaic(self, cutouts):
    
        td = tempfile.mkdtemp()
        input_dir = '{directory}/input'.format(directory=td)
        output_dir = '{directory}/output'.format(directory=td)
        self.make_dir(input_dir)
    
        try:
            for i, c in enumerate(cutouts):
    
                with open('{directory}/{name}.fits'.format(directory=input_dir, name=i), 'wb') as tmp:
                    tmp.write(bytes(c))
            os.listdir(td)
            montage_wrapper.mosaic(input_dir, output_dir)
    
            with open('{outdir}/mosaic.fits'.format(outdir=output_dir), 'rb') as f:
    
                merged = f.read()
    
        finally:
            shutil.rmtree(output_dir)
            shutil.rmtree(input_dir)
            shutil.rmtree(td)
    
        return merged
