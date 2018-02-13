from ftplib import FTP

ncbi_ftp = 'ftp.ncbi.nlm.nih.gov'
orange = 'genomes/Citrus_sinensis'
fasta = 'csi_ref_Csi_valencia_1.0_chr1.fa.gz'

ftp = FTP(ncbi_ftp, user='anonymous', passwd='gribskov@purdue.edu')
ftp.cwd(orange)  # change working directory

for file in ftp.nlst():
    if file.startswith('CHR'):
        print('\nstarting', file)
        ftp.cwd(file)
        for datafile in ftp.nlst():
            print(datafile, end=' ')
            if datafile.endswith('.fa.gz'):
                print('--> retrieving {} ... '.format(datafile), end=' ')
                try:
                    ftp.retrbinary("RETR " + datafile, open(datafile, 'wb').write)
                    print('done')
                except:
                    print('Error retrieving {}'.format(fasta))
            else:
                print('   skip')
        ftp.cwd('..')

exit(0)
