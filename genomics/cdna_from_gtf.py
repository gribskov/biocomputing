"""=================================================================================================
get cDNA sequences from genome fa file based on Ensembl GTF

Michael Gribskov     03 July 2024
================================================================================================="""
import sys

sys.path.append('../gff')
from gtf2gff import Gtf

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    gtf_file = 'A:/mrg\Dropbox/21dog_dlbcl/reference data/cfam112/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf'
    gtf = Gtf(gtf_file)

    while gtf.next():
        gtf.parse()
        # gtf.attribute2gff()
        gtf.add_ensemble_attributes()
        print(f'{gtf.parsed}\n')

    exit(0)
