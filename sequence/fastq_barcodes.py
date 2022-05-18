"""=================================================================================================
get the barcodes from a set of fastq files.  for the newer illumina double barcodes with ID lines
like
@M03613:204:000000000-JLYF9:1:1119:18149:1493 1:N:0:TCGCGCAT+TGAATTGG

Michael Gribskov     18 May 2022
================================================================================================="""
import glob
import sys
import os.path
from fastq import Fastq

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    suffix = '.fastq'

    target = sys.argv[1]
    if not target.endswith('/'):
        target += '/'
    target += f'*{suffix}'

    codes = {}
    counts = {}
    for fastqfile in glob.glob(target):
        print(f'{fastqfile}')
        fastq = Fastq(fastqfile)

        # get the sample name from the filename
        sample = os.path.basename(fastqfile).replace(suffix, '')
        samplefields = sample.split('_')
        sample = '_'.join(samplefields[:2])

        while fastq.next():
            spot, codestring = fastq.id.split()

            barcodefields = codestring.split(':')
            barcode = barcodefields[3].replace('+', '')
            print(f'{barcode}', end='\t')
            print(f'{fastq.id}\t{sample}')

            if sample in codes:
                if barcode not in codes[sample]:
                    codes[sample].append(barcode)
                    counts[sample][barcode] = 1
                else:
                    counts[sample][barcode] += 1
            else:
                codes[sample] = [barcode]
                counts[sample] = {barcode:1}

    for sample in codes:
        print(f'{sample}')
        for barcode in codes[sample]:
            print(f'\t{barcode}\t{counts[sample][barcode]}')

    exit(0)
