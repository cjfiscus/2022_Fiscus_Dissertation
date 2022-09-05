# combine files
#
# cjfiscus

base_name = "../../results/feat_abund/parts2/part_"

content = []
LineNo = 0
start = 2 # first col w/ abundance estimates in kmer table
end = 1048 # last col of kmer table

# collect data into strings with each line as item in content[]
for i in range(start,(end + 1)):
    FILENAME = str(base_name) + str(i) + ".txt"
    File = open(FILENAME, "r")
    LineNo = 0

    for Line in File:
        Value = Line.split("\t")[1]

        if i == start:
            content.append(Value)

        else:
            content[LineNo] = str(content[LineNo]) + "\t" + str(Value)
            LineNo = LineNo + 1

    File.close()

# print out combined data
for i in range(len(content)):
    print(str(content[i]))
