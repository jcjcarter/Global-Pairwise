import random, numpy, collections, comp182

def global_pairwise_alignment(sequen1, sequen2, m):
    """Finds the best alignment sequence for two inputted
    strings.

    Arguments:
    sequen1 -- string
    sequen2 -- string
    m -- scoring matrix

    Returns:
    Best alignment sequence for the two inputted strings.
    """
    i = len(sequen1) 
    j = len(sequen2)  
    possmoves = [(-1, -1), (-1, 0), (0, -1)] #List of possibly moves in the matrix for trace back
    xstr = ''
    ystr = ''
    score = numpy.zeros((i+1,j+1)) #initializes the matrix
    score[0][0] = 0
    counter = 1
    
    for a in range(1, i+1): #scores the rows 
        
        score[a][0] = max([score[a-1][0] + m[sequen1[counter-1]]['-']]) 
        counter += 1
    counter = 1
    for a in range(1, j+1): #scores the column
        
        score[0][a] = max([score[0][a-1] + m['-'][sequen2[counter-1]]]) 
        counter += 1
    r = 0
    for rowcounter in range(1, i+1): #gets the diagonal, column, and row
        c = 0
        for colcounter in range(1,j +1):
            across = score[rowcounter-1][colcounter-1] + m[sequen1[rowcounter-1]][sequen2[colcounter-1]]
            above   = score[rowcounter-1][colcounter] + m[sequen1[rowcounter-1]]['-']
            left = score[rowcounter][colcounter-1] + m['-'][sequen2[colcounter-1]]
            score[rowcounter][colcounter] = max(across, above, left) 
            c += 1
        r += 1
    while (i > 0) and (j > 0): #finds the highest value in the scoring matrix
        across = score[i-1][j-1] + m[sequen1[i-1]][sequen2[j-1]]
        above = score[i-1][j] + m[sequen1[i-1]]['-']
        left = score[i][j-1] + m['-'][sequen2[j-1]]
        neigh = [across, above, left]
        maxneigh = max(neigh)

        find = neigh.index(maxneigh)
        if possmoves[find][0] == -1:
            i -= 1
            xstr = sequen1[i] + xstr
        else:
            xstr = '-' + xstr
        if possmoves[find][1] == -1:
            j -= 1
            ystr = sequen2[j] + ystr
        else:
            ystr = '-' + ystr
    while i > 0: #adds the remaining letters left in the sequence
        i -= 1
        xstr = sequen1[i] + xstr
        ystr = '-' + ystr
    while j > 0:
        j -= 1
        ystr = sequen2[j] + ystr
        xstr = '-' + xstr

        
    return xstr, ystr

def scoreMatrix():
    """Constructs a score matrix.

    Arguments:
    None

    Returns:
    The ACTG- score matrix"""

    codones = 'ACTG-'
    score = {} #the score matrix that will be returned
    for i in codones:
        score[i] = {}
        for j in codones:
            value = 2
            if j == i:
                value = 6
            if (i == '-') or (j == '-'):
                value = -4
            score[i][j] = value
    return score

def test_global_pairwise_alignment():
    """Tests the functionality of the function
    global_pairwise_alignment.

    Arguments:
    None

    Returns:
    Results of the test."""

    m = scoreMatrix()
    string1 = 'ACTGGG-'
    string2 = 'ACT-G-'
    return global_pairwise_alignment(string1, string2, m)
    

def local_pairwise_alignment(sequen1, sequen2, m):
    """Implements...

    Arguments:
    None

    Returns:
    Substrings of X and Y whose global alignment (as defined
    by the matrix) is maximal among all global alignments
    of all substrings of X and Y."""
    
    i = len(sequen1) 
    j = len(sequen2)  
    possmoves = [(-1, -1), (-1, 0), (0, -1)] #possible moves in the matrix
    xstr = ''
    ystr = ''
    score = numpy.zeros((i+1,j+1))
    
    for a in range(1, i+1): #scores the rows    
        
        score[a][0] = findMax([score[a-1][0] + m[sequen1[a-1]]['-']]) 
        
    for a in range(1, j+1): #scores the columns
        score[0][a] = findMax([score[0][a-1] + m['-'][sequen2[a-1]]]) 
        
    r = 0
    for row in range(1, i+1): # gets score for the diagonal and remaining columns/rows
        c = 0
        for col in range(1,j +1):
            across = score[row-1][col-1] + m[sequen1[row-1]][sequen2[col-1]]
            above   = score[row-1][col] + m[sequen1[row-1]]['-']
            left = score[row][col-1] + m['-'][sequen2[col-1]]
            score[row][col] = findMax([across, above, left]) 
            
    xx = -1 #starts at the end of the sequence
    xy = -1
    ml = None
    for i in range(i+1):
        for j in range(j+1):
            if (ml == None) or (score[i, j] > ml):
                ml = score[i, j]
                xx = i
                xy = j
    i = xx #finds the highest value for the string
    j = xy
    
    while (i > 0) and (j > 0): #performs the backtrace to find the highest score
        across = score[i-1][j-1] + m[sequen1[i-1]][sequen2[j-1]]
        above = score[i-1][j] + m[sequen1[i-1]]['-']
        left = score[i][j-1] + m['-'][sequen2[j-1]]
        neigh = [across, above, left]
        maxneigh = max(neigh)
        if score[i, j] != maxneigh:
            ### Local alignment, S[i, j] should be 0
            if score[i, j] != 0:
                print "Problem at:", score[i,j]
            return xstr, ystr
        find = neigh.index(maxneigh)
        if possmoves[find][0] == -1:
            i -= 1
            xstr = sequen1[i] + xstr
        else:
            xstr = '-' + xstr
        if possmoves[find][1] == -1:
            j -= 1
            ystr = sequen2[j] + ystr
        else:
            ystr = '-' + ystr
    
    while i > 0:
        i -= 1
        xstr = sequen1[i] + xstr
        ystr = '-' + ystr
    while j > 0:
        j -= 1
        ystr = sequen2[j] + ystr
        xstr = '-' + xstr

    return xstr, ystr

def findMax(listvalues):
    """Finds the  highest value in the list given.

    Arguments:
    listvalues -- list of values

    Returns:
    0 if highest value is less than 0  or the highest found value."""
    fmax = max(listvalues)
    if fmax < 0:
        fmax = 0
    return fmax
def test_local_pairwise_alignment():
    """Tests the functionality of the function
    local_pair_alignment.

    Arguments:
    None

    Returns:
    Results of the test."""
    m = scoreMatrix()
    string1 = 'ACTGGG-'
    string2 = 'ACT-G-'
    return local_pairwise_alignment(string1, string2, m)

def generate_null_distribution(x,y,M,t,numiter):
    """Functions takes the two sequences x and y, and repeats the
    following process for numiter times: generate a random permutation
    y' of the sequence y, align x with y' using the score matrix M,
    and compute the score of the best alignment. The function should
    plot the distribution of the numiter scores.

    Arguments:
    x -- sequence
    y -- sequence
    M -- score matrix
    t -- an integer, 1 or 0, 0 for local and 1 for global
    numiter -- an integer

    Returns:
    The mean of the numiter scores, the standard deviation of the numiter scores,
    the value z, and plot distribution."""

    listscored = []
    if t == 0: #from the input determines what input to use
        fucn = local_pairwise_alignment
    else:
        fucn = global_pairwise_alignment
    xfucn, yfucn = fucn(x,y,M)
    
    se = score_alignment(xfucn, yfucn, M) #produces the actual score
    print ''
    
    for a in range(numiter): # generates random strings and scores them
        stringY = permute_string(y)
        xstr, ystr = fucn(x, stringY, M)
        result = score_alignment(xstr, ystr, M)
        
        listscored.append(result)
    print 'Scored', listscored
    cd = convertDic(listscored)
    scaleDict(cd)
    comp182.plot_dist_linear(cd, 'Distribution', 'Scores', 'Fraction',filename='Null Distribution' )
    av = numpy.average(listscored) #computes the mean
     
    sd = numpy.std(listscored) #computes standard deviation
    
    #comp182.show()
    return 'Average:', numpy.average(listscored),'Standard Deviation:', sd,'Z-Value:', (se - av)/sd

def scaleDict(d):
    """Scales the dictionary values by dividing all entries by
    the sum of the values.

    Arguments:
    d -- dictionary to be normalized

    Returns:
    None
    """
    s = float(sum(d.values()))
    if s != 0.0:
        for k in d:
            d[k] = float(d[k]) / s

def convertDic(l):
    """
    Convert a list of values to a distribution.

    Arguments:
    l -- list

    Returns:
    Dictionary with unique elements of l as keys and the number of
    times that key occurs in l as values
    """
    d = collections.defaultdict(int)
    for i in l:
        d[i] += 1
    return d
    

def test_generate_null_distribution():
    """Tests the functionality of the function generate_null_distribution.

    Arguments:
    None

    Returns:
    Results of the test."""
    M = scoreMatrix()
    sx = 'ACTG--ACTCA'
    sy = 'AGTC--ATGC'
    print generate_null_distribution(sx,sy,M,1,20)
    print ''
    print generate_null_distribution(sx,sy,M,0,20)

def permute_string(s):
    """
    Return a new string with the characters in s randomly permuted.

    Arguments:
    s -- string

    Returns:
    Random permutation of s
    """
    charlist = list(s)
    random.shuffle(charlist)
    newstr = "".join(charlist)
    return newstr

def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    with open(filename) as f:
        p = f.read()
    p = p.rstrip()
    return p


def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    M = {}
    with open(filename) as f:
        ykeys = f.readline()
        ykeychars = ykeys.split()
        for line in f.readlines():
            vals = line.split()
            xkey = vals.pop(0)
            M[xkey] = {}
            for ykey, val in zip(ykeychars, vals):
                M[xkey][ykey] = int(val)
    return M

def runEyeless():
    """Run the file eyeless through the functions global_pairwise_alignment
    and local_pairwise_alignment.

    Argument:
    None

    Returns:
    Results from the run of two functions."""
    sm = read_scoring_matrix('PAM50.txt')
    hp = read_protein('HumanEyelessProtein.txt')
    fy = read_protein('FruitflyEyelessProtein.txt')
    print "Local Alignment Results..."
    sx, sy = local_pairwise_alignment(hp, fy, sm)
    print ''
    print "Human", sx
    print ''
    print "Fly", sy
    print ''

    print "Global Alignment Results..."
    sx, sy = global_pairwise_alignment(hp, fy, sm)
    print ''
    print "Human", sx
    print ''
    print "Fly", sy

def distributionEyeless(run):
    """Run the file eyeless through the functions global_pairwise_alignment
    and local_pairwise_alignment and produces the null distribution.

    Argument:
    run -- number of times the loop is executed in distribution.

    Returns:
    Results from the run of two functions."""
    
    sm = read_scoring_matrix('PAM50.txt')
    hp = read_protein('HumanEyelessProtein.txt')
    fy = read_protein('FruitflyEyelessProtein.txt')
    print "Local Alignment Results..."
    print generate_null_distribution(hp, fy, sm, 0, run)
    print ''

    print "Global Alignment Results..."
    print generate_null_distribution(hp, fy, sm, 1, run)
    

def score_alignment(xstr, ystr, M):
    """Produces the alignment's score based off of the scoring matrix.

    Arguments:
    xstr -- sequence
    ystr -- sequence
    M -- scoring matrix

    Returns:
    Score for the alignment
    """
    se = 0
    for a, b in zip(xstr, ystr):
        se += M[a][b]
    return se
