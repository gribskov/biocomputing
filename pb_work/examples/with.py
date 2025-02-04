"""=================================================================================================
example of reading file using with

Michael Gribskov     30 January 2025
================================================================================================="""

# --------------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # read whole file (watch out for memory issues)
    with open("wrapper.py") as my_file:
        print(my_file.read())

    # read line by line, not the extra indentation
    with open("wrapper.py") as my_file:
        line = my_file.readline()
        while line:
            print(line.rstrip())
            line = my_file.readline()

    # my_file.close() is not needed when using with

    after()

    exit(0)
