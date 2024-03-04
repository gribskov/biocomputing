"""=================================================================================================
query phobius server at http://phobius.sbc.su.se/

Michael Gribskov     18  2018
================================================================================================="""
import requests
from bs4 import BeautifulSoup

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    phobius = 'https://phobius.sbc.su.se/'
    submit = 'cgi-bin/predict.pl'

    bri1 = '''>AEE87069.1 Leucine-rich receptor-like protein kinase family protein [Arabidopsis thaliana]
MKTFSSFFLSVTTLFFFSFFSLSFQASPSQSLYREIHQLISFKDVLPDKNLLPDWSSNKNPCTFDGVTCR
DDKVTSIDLSSKPLNVGFSAVSSSLLSLTGLESLFLSNSHINGSVSGFKCSASLTSLDLSRNSLSGPVTT
LTSLGSCSGLKFLNVSSNTLDFPGKVSGGLKLNSLEVLDLSANSISGANVVGWVLSDGCGELKHLAISGN
KISGDVDVSRCVNLEFLDVSSNNFSTGIPFLGDCSALQHLDISGNKLSGDFSRAISTCTELKLLNISSNQ
FVGPIPPLPLKSLQYLSLAENKFTGEIPDFLSGACDTLTGLDLSGNHFYGAVPPFFGSCSLLESLALSSN
NFSGELPMDTLLKMRGLKVLDLSFNEFSGELPESLTNLSASLLTLDLSSNNFSGPILPNLCQNPKNTLQE
LYLQNNGFTGKIPPTLSNCSELVSLHLSFNYLSGTIPSSLGSLSKLRDLKLWLNMLEGEIPQELMYVKTL
ETLILDFNDLTGEIPSGLSNCTNLNWISLSNNRLTGEIPKWIGRLENLAILKLSNNSFSGNIPAELGDCR
SLIWLDLNTNLFNGTIPAAMFKQSGKIAANFIAGKRYVYIKNDGMKKECHGAGNLLEFQGIRSEQLNRLS
TRNPCNITSRVYGGHTSPTFDNNGSMMFLDMSYNMLSGYIPKEIGSMPYLFILNLGHNDISGSIPDEVGD
LRGLNILDLSSNKLDGRIPQAMSALTMLTEIDLSNNNLSGPIPEMGQFETFPPAKFLNNPGLCGYPLPRC
DPSNADGYAHHQRSHGRRPASLAGSVAMGLLFSFVCIFGLILVGREMRKRRRKKEAELEMYAEGHGNSGD
RTANNTNWKLTGVKEALSINLAAFEKPLRKLTFADLLQATNGFHNDSLIGSGGFGDVYKAILKDGSAVAI
KKLIHVSGQGDREFMAEMETIGKIKHRNLVPLLGYCKVGDERLLVYEFMKYGSLEDVLHDPKKAGVKLNW
STRRKIAIGSARGLAFLHHNCSPHIIHRDMKSSNVLLDENLEARVSDFGMARLMSAMDTHLSVSTLAGTP
GYVPPEYYQSFRCSTKGDVYSYGVVLLELLTGKRPTDSPDFGDNNLVGWVKQHAKLRISDVFDPELMKED
PALEIELLQHLKVAVACLDDRAWRRPTMVQVMAMFKEIQAGSGIDSQSTIRSIEDGGFSTIEMVDMSIKE
VPEGKL'''

    command = {'protseq': bri1, 'format': 'plp'}
    print(f'Sending query to {phobius + submit}')
    response = requests.post(phobius + submit, data=command)

    # find the image tag with beautiful soup, the image shows the prediction probabilities
    soup = BeautifulSoup(response.content, 'html.parser')
    img = soup.find('img')

    # plot address is in the image tag source, beginning at character 3
    plot = requests.get(phobius + img['src'][3:], stream=True)
    print(f"status:{plot} prediction plot: {img['src']}")
    with open('phobius.png', 'wb') as image:
        # write out plot in binary format
        image.write(plot.content)

    exit(0)
