from ftplib import FTP

# https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Citrus_sinensis/latest_assembly_versions/GCA_000695605.1_Citrus_sinensis_v1.0/
ncbi_ftp = 'ftp.ncbi.nlm.nih.gov'
genome_plant = 'genomes/genbank/plant'
orange = 'Citrus_sinensis/latest_assembly_versions/GCA_000695605.1_Citrus_sinensis_v1.0'

ftp = FTP(ncbi_ftp, user='anonymous', passwd='gribskov@purdue.edu')
cisen = '{}/{}'.format(genome_plant, orange)
ftp.cwd(cisen)  # change working directory

# directory listing
filelist = []
ftp.dir(filelist.append)
print(*filelist, sep='\n')

print('\nTransferring files')
for file in sorted(ftp.nlst()):
    print('\t{}'.format(file), end=' ')
    if file.startswith('GCA') and file.endswith('.gz'):
        print('--> retrieving', end=' ')
        try:
            ftp.retrbinary("RETR " + file, open(file, 'wb').write)
            print('--> done')
        except:
            print('--> error')
    else:
        print('--> skip')

exit(0)
