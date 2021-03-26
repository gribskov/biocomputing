"""=================================================================================================
Examples of subprocess

Michael Gribskov     16 March 2021
================================================================================================="""
import subprocess as sub

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':

   # run OS command, output to stdout
    sub.run(['dir'], shell=True)


   # run OS command with options, output to stdout
   #  sub.run(['dir', 'blast.py'], shell=True)
   #  sub.run(['dir', '/w'], shell=True)

   # OS command, with result captured as text
   #  result = sub.run(['dir', '/w'], shell=True, capture_output = True, text = True )
   #  print("stdout:", result.stdout)
   #  print("stderr:", result.stderr)

   # OS command, with result captured as text, with error
   #  result = sub.run(['dir', '/ww'], shell=True, capture_output = True, text = True)
   #  print("stdout:", result.stdout)
   #  print("stderr:", result.stderr)

   # check if command was successful (returned status=0)
   #  try:
   #      result = sub.run(['dir', '/ww'], shell=True, check=True)
   #  except sub.CalledProcessError:
   #      print( 'Command failed')

    first = sub.run(['type', 'blast.py'], stdout=sub.PIPE, text=True, shell=True)
    output = sub.check_output(['find', '"def"'], text=True, stdin=first.stdout, shell=True)
    print('x', output.stdout)
    # output = sub.Popen(('find', '"."'), stdin=first.stdout)
    # print(output.stdout)

    # dir = sub.run(['dir'], stdout=sub.PIPE, shell=True)
    # sub.run(['type'], stdin=dir.stdout, universal_newlines=True, shell=True)
    # sub.run(['find', '10/26/2018'], stdin=dir.stdout, universal_newlines=True, shell=True)
    # print(date.stdout)


    # output = subprocess.check_output(('grep', 'process_name'), stdin=ps.stdout)


    # keep internal input, convert byte result to string, windows specific newline split
    # result = sub.run(['dir', 'blast*'], stdout=sub.PIPE, shell=True)
    # for line in result.stdout.decode().split('\r\n'):
    #     print(line)

    # # automatic decoding using universal_newlines=True
    # result = sub.run(['dir', ''], stdout=sub.PIPE, universal_newlines=True, shell=True)
    # for line in result.stdout.split('\r\n'):
    #     print(line)
    #
    # # keep internal input, convert byte result to string, windows specific newline split
    # # with error in standard error
    # result = sub.run(['dir', '-f'], stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
    # print('stderr:', result.stderr.decode())

    exit(0)
