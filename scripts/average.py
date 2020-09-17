import sys
import os.path

if( len(sys.argv) < 2):
    print("Error. Usage: python3 average.py [data file]")
    sys.exit(0)
try:
    filename = sys.argv[1]

    #open and read file
    INFILE = open(filename,"r")
    data = INFILE.readlines()
    if(len(data) < 1):
        print("No data.")
        sys.exit(0)
    INFILE.close()

    #remove trailing newline entry
    if(data[-1] == ""):
        del data[-1]

    #calculate average
    sum = 0
    count = 0
    for i in data:
        sum += float(i.rstrip())
        count += 1
    ave = sum/count

    #generate output file name
    output_path = filename.split("/")
    output_filename = output_path[-1].rstrip().split(".")[0] + "_ave.txt"
    output_path = output_path[:-1]
    output_path.append(output_filename)
    output_file = "/".join(output_path)

    #write to file
    OUTFILE = open(output_file,"w")
    OUTFILE.write(str(ave))
    OUTFILE.close()

except FileNotFoundError:
    print("File {} not found".format(filename))
except ValueError:
    print("Something wrong with the values in the file.")

