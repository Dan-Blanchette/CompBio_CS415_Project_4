# from IPython.display import HTML
# shell = get_ipython()

# def adjust_font_size():
#     display(HTML('''<style>
#         body{
#             font-size: 32px;
#         }
#     '''))

# if adjust_font_size not in shell.events.callbacks['pre_execute']:
#     shell.events.register('pre_execute', adjust_font_size)

def printmatrix(m,pad = 4): 
    for r in m:
        for d in r:
            print(f"{str(d):>{pad}}", end = " ")
        print()

# Dictionary to map nucleic acid codes to matrix indices
aminoDictionary = {'a':0, 'r':1, 'n':2, 'd':3, 'c':4, 'q':5, 'e':6, 'g':7,'h':8, 
'i':9, 'l':10, 'k':11, 'm':12, 'f':13, 'p':14, 's':15,'t':16, 'w':17, 'y':18, 
'v':19}
#Gap penalty
gappenalty = -8
#blosum[0][0] = 5
#blosum50[0][1] = -2
#blosum50[0][2] = -1
#... very slow
#blosum50 matrix
blosum50 = [[]]
blosum50[0] = [5,-2,-1,-2,-1,-1,-1,0,-2,-1,-2,-1,-1,-3,-1,1,0,-3,-2,0]
blosum50.append([-2,7,-1,-2,-4,1,0,-3,0,-4,-3,3,-2,-3,-3,-1,-1,-3,-1,-3])
blosum50.append([-1,-1,7,2,-2,0,0,0,1,-3,-4,0,-2,-4,-2,1,0,-4,-2,-3])
blosum50.append([-2,-2,2,8,-4,0,2,-1,-1,-4,-4,-1,-4,-5,-1,0,-1,-5,-3,-4])
blosum50.append([-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1])
blosum50.append([-1,1,0,0,-3,7,2,-2,1,-3,-2,2,0,-4,-1,0,-1,-1,-1,-3])
blosum50.append([-1,0,0,2,-3,2,6,-3,0,-4,-3,1,-2,-3,-1,-1,-1,-3,-2,-3])
blosum50.append([0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4])
blosum50.append([-2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4])
blosum50.append([-1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4])
blosum50.append([-2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1])
blosum50.append([-1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3])
blosum50.append([-1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1])
blosum50.append([-3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1])
blosum50.append([-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3])
blosum50.append([1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2])
blosum50.append([0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0])
blosum50.append([-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,2,-3])
blosum50.append([-2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1])
blosum50.append([0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5])

seq1 = "heagawghee"
seq2 = "pawheae"
#gappenalty = -8
scoringmatrix = [[00 for i in range(0,len(seq1)+1)] for j in range(0,len(seq2)+1)]
directionmatrix = [["." for i in range(0,len(seq1)+1)] for j in range(0,len(seq2)+1)]
# printmatrix(scoringmatrix)

for i in range(0,len(scoringmatrix[0])):
    scoringmatrix[0][i] =  i*gappenalty
    directionmatrix[0][i] = '\u2190'
for i in range(0,len(scoringmatrix)):
    scoringmatrix[i][0] = i*gappenalty
    directionmatrix[i][0]  = '\u2191'

for r in range(1,len(scoringmatrix)):
    for c in range(1,len(scoringmatrix[r])):
        vert = scoringmatrix[r-1][c] + gappenalty
        horz = scoringmatrix[r][c-1] + gappenalty
        diag = scoringmatrix[r-1][c-1]
        #print(seq1[r-1],seq2[c-1])
# seq1 = "heagawghee"
#seq2 = "pawheae"
        char1 = seq1[c-1]
        char2 = seq2[r-1]
        index1 = aminoDictionary[char1]
        index2 = aminoDictionary[char2]
        diag += blosum50[index2][index1]
      #  diag += blosum50[aminoDictionary[seq2[r-1]]][aminoDictionary[seq1[c-1]]]
       # if(seq1[c-1] == seq2[r-1]):
       #     diag += 1
       # else:
       #     diag -= 1
        scoringmatrix[r][c] = max(vert,horz,diag)
        #directionmatrix[r][c] = "v or h or d"
        if diag >= horz and diag >= vert:
            directionmatrix[r][c] = "\u2196"
        if horz > diag and horz > vert:
            directionmatrix[r][c] = "\u2190"
        if vert > diag and vert > horz:
            directionmatrix[r][c] = "\u2191"

printmatrix(scoringmatrix,4)
print(seq1)
print(seq2)

printmatrix(directionmatrix,0)

aligned1 = ""
aligned2 = ""
pos1 = len(seq1)
pos2 = len(seq2)

while pos1 > 0 or pos2 > 0:
    if(directionmatrix[pos2][pos1] == '\u2196'): # diag
        aligned1 += (seq1[pos1-1])
        aligned2 += (seq2[pos2-1])
        pos1 -=1 
        pos2 -=1
    else:
        if(directionmatrix[pos2][pos1] == '\u2190'): # horz
            aligned1 += (seq1[pos1-1])
            aligned2 += ('-')
            pos1 -=1
        else: 
            if(directionmatrix[pos2][pos1] == '\u2191'): # vert
                aligned1 += ('-')
                aligned2 += (seq2[pos2-1]) 
                pos2 -=1
            else: 
                print("Something is wrong!")
aligned1 = aligned1[::-1]
aligned2 = aligned2[::-1]
print(aligned1)
print(aligned2)