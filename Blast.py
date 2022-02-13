""""
Team member:
nouran ahmed ibrahim
sama medhat farouk
menna tullah abdelhalem
nada yoseef abdelmo3z
yomnna atef abdelmon3m
abeer husien mohamed
"""
from tqdm import tqdm

query_sequence = open("query", 'r').readline().upper()
seeds = []
database_seqs = {}
amino = {}
results = {}
ranking = []
blosum62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}

similar_amino_acids = {
    # 9
    'C': ['C'],
    # 8
    'H': ['H'],
    # 11
    'W': ['W'],
    # 7
    'Y': ['Y', 'P'],
    'P': ['P', 'Y'],
    # 5
    'T': ['T', 'E', 'Q', 'R', 'K', 'M'],
    'E': ['E', 'T', 'Q', 'R', 'K', 'M'],
    'Q': ['Q', 'E', 'T', 'R', 'K', 'M'],
    'R': ['R', 'K', 'M', 'Q', 'E', 'T'],
    'K': ['K', 'M', 'R', 'Q', 'E', 'T'],
    'M': ['M', 'K', 'R', 'Q', 'E', 'T'],
    # 4
    'A': ['A', 'I', 'L', 'S', 'V'],
    'I': ['I', 'A', 'L', 'S', 'V'],
    'L': ['L', 'A', 'I', 'S', 'V'],
    'S': ['S', 'A', 'L', 'I', 'V'],
    'V': ['V', 'A', 'I', 'L', 'S'],
    # 6
    'G': ['G', 'N', 'D', 'F'],
    'N': ['N', 'G', 'D', 'F'],
    'D': ['D', 'G', 'N', 'F'],
    'F': ['F', 'D', 'N', 'G'],
}


"""
remove_low_complexity_region:
Remove low complexity region by going through the protein and looking for a similar triple sequence.
Parameters: protein query sequence
Return: -
"""
def remove_low_complexity_region(protein_seq=query_sequence):
    iterator = 0  # iterate on protein seq

    while iterator < len(protein_seq) - 5:
        start_index = -1
        end_index = -1

        for i in range(iterator, len(protein_seq) - 5):
            current_amino_acid = protein_seq[i]
            next_amino_acid = protein_seq[i + 1]

            if (protein_seq[i + 2] in similar_amino_acids[current_amino_acid] and
                    protein_seq[i + 4] in similar_amino_acids[current_amino_acid] and
                    protein_seq[i + 3] in similar_amino_acids[next_amino_acid] and
                    protein_seq[i + 5] in similar_amino_acids[next_amino_acid]):

                if start_index == -1:  # first three repeated region
                    start_index = i

                end_index = i + 5

            else:
                if start_index != -1:
                    break

        if start_index != -1:  # if repeated region is found in sequence
            number = end_index - start_index + 1
            protein_seq = protein_seq.replace(protein_seq[start_index: end_index + 1], 'X' * number)
            iterator = end_index + 1

        else:  # if no repeated region in sequence
            iterator = len(protein_seq) - 5

    return protein_seq


"""
read_database:
Read the database sequences from file and append them in dictionary (key is id and value is sequence).
Parameters: -
Return: -
"""
def read_database():
    file = open("database", 'r')
    lines = file.readlines()
    file.close()

    for line_idx in range(0, len(lines), 2):
        id = lines[line_idx].split('\n')[0]
        seq = lines[line_idx + 1].split('\n')[0].upper()
        database_seqs[id] = seq

"""
create_initial_seeds:
Make word with size w and append it in list if it's score is greater than or equal word threshold
Parameters: -
Return: -
"""
def create_initial_seeds():
    for i in range(len(query_sequence) - w + 1):
        end = i + w - 1
        word = query_sequence[i:end + 1]
        if query_sequence[i] != 'X' and query_sequence[i + 1] != 'X' and query_sequence[i + 2] != 'X':
            score = calc_score(word, word)
            if score >= word_Threshold:
                seeds.append([word, [i, end, score]])


"""
calc_score:
Calculate the score of two string 
Parameters: string, string
Return: Score
"""
def calc_score(str1, str2):
    score = 0
    for char_idx in range(len(str1)):
        if (str1[char_idx], str2[char_idx]) in blosum62.keys():
            score += blosum62[(str1[char_idx], str2[char_idx])]
        else:
            score += blosum62[(str2[char_idx], str1[char_idx])]

    return score


"""
make_permutation:
Find all the permutations for each word in the seeds list and append it in the seeds list if its value is greater than or equal to word threshold.
Parameters: -
Return: -
"""
def make_permutation():
    find_last_word = False
    last_word = seeds[len(seeds) - 1]
    for word in seeds:
        if find_last_word:
            break
        for word_char_index in range(len(word[0])):
            if word == last_word:
                find_last_word = True
            original_word_char = word[0][word_char_index]
            old_amino_score = calc_score(original_word_char, original_word_char)
            for amino in similar_amino_acids:
                if amino == original_word_char:
                    continue
                new_amino_score = calc_score(amino, original_word_char)
                score = word[1][2]
                final_score = score + new_amino_score - old_amino_score
                if final_score >= word_Threshold:
                    edit_word = change_amino(word[0], amino, word_char_index)
                    start = word[1][0]
                    end = word[1][1]
                    seeds.append([edit_word, [start, end, final_score]])

"""
change_amino:
Replace the letter located in the index with the letter to be replaced.
Parameters: original word and its index, character to replace.
Return: new word.
"""
def change_amino(original_word, new_char, position):
    if position == 0:
        return new_char + original_word[1:]
    elif position == (len(original_word) - 1):
        return original_word[0:len(original_word) - 1] + new_char
    else:
        return original_word[0:position] + new_char + original_word[position + 1:]


"""
HSP:
calculate the high score segment pair and return score index and max score .
Parameters: id of the database sequence, start and end index of database sequence, start and end index of query sequence, fist iteration (Boolean), max score.
Return: max score, score index.
"""
def HSP(id, Start_DB, End_DB, Start_Seq, End_Seq, first_iteration, maxScore, ScoreIdx):
    DBSeq = database_seqs[id]
    Hit_Seq = query_sequence[Start_Seq:End_Seq + 1]
    Hit_DB = DBSeq[Start_DB:End_DB + 1]
    Score = calc_score(Hit_Seq, Hit_DB)
    if first_iteration:
        maxScore = Score
        ScoreIdx = [Start_DB, End_DB, Start_Seq, End_Seq]

    if Start_DB - 1 >= 0 and End_DB + 1 <= len(DBSeq) - 1 and Start_Seq - 1 >= 0 and End_Seq + 1 <= len(query_sequence) - 1:
        if maxScore - Score <= 5:
            if maxScore < Score:
                maxScore = Score
            ScoreIdx = [Start_DB, End_DB, Start_Seq, End_Seq]
            return HSP(id, Start_DB - 1, End_DB + 1, Start_Seq - 1, End_Seq + 1, False, maxScore, ScoreIdx)
        else:
            return maxScore, ScoreIdx
    elif (Start_DB - 1 < 0 or Start_Seq - 1 < 0) and End_DB + 1 <= len(DBSeq) - 1 and End_Seq + 1 <= len(
            query_sequence) - 1:
        if maxScore - Score <= 5:
            if maxScore < Score:
                maxScore = Score
            ScoreIdx = [Start_DB, End_DB, Start_Seq, End_Seq]
            return HSP(id, Start_DB, End_DB + 1, Start_Seq, End_Seq + 1, False, maxScore, ScoreIdx)
        else:
            return maxScore, ScoreIdx

    elif Start_DB - 1 >= 0 and Start_Seq - 1 >= 0 and (
            End_DB + 1 > len(DBSeq) or End_Seq + 1 > len(query_sequence)):
        if maxScore - Score <= 5:
            if maxScore < Score:
                maxScore = Score
            ScoreIdx = [Start_DB, End_DB, Start_Seq, End_Seq]
            return HSP(id, Start_DB - 1, End_DB, Start_Seq - 1, End_Seq, False, maxScore, ScoreIdx)
        else:
            return maxScore, ScoreIdx
    else:
        return maxScore, ScoreIdx


"""
Find_Similar_seqs:
Find all Hits in each sequence that have score greater than HSP_Thershold.
Parameters: -
Return: -
"""
def find_Similar_seqs():
    for id in database_seqs:
        DBSeq = database_seqs[id]
        total_score = 0
        seq_hits = []
        for i in range(0, len(DBSeq) - w + 1):
            word = DBSeq[i: i + 3]
            for seed in seeds:
                if seed[0] == word:
                    hit_score, idx = HSP(id, i, i + 2, seed[1][0], seed[1][1], True, 0, [])
                    star_DBHit_Idx = idx[0]
                    end_DBHit_Idx = idx[1]
                    star_QHit_Idx = idx[2]
                    end_QHit_Idx = idx[3]
                    if hit_score >= HSP_Thershold:
                        total_score += hit_score
                        hit = [star_DBHit_Idx, end_DBHit_Idx, star_QHit_Idx, end_QHit_Idx, hit_score]
                        seq_hits.append(hit)
        if len(seq_hits) != 0:
            results[id] = [seq_hits, total_score]
            ranking.append([total_score, id])
    ranking.sort(reverse=True)


"""
print_Similar_Seqs:
Display Hits in each sequence that have score greater than HSP_Thershold.
Sequences are ordered descending based on its total score.
It display hit sequence in database sequence and in query sequence and its score.
Parameters: -
Return: -
"""
def print_Similar_Seqs():
    for seqR in ranking:
        seqScore = seqR[0]
        seqId = seqR[1]
        print("\nSeq ", seqId, " has ", len(results[seqId][0]), " hits with total score", seqScore)
        print('Hits:')
        hit_num = 1
        for hit in results[seqId][0]:
            DBHit = database_seqs[seqId][hit[0]:hit[1] + 1]
            QHit = query_sequence[hit[2]:hit[3] + 1]
            print("Hit ", hit_num)
            print("Hit seq in database starts from ", hit[0], " to ", hit[1], " and hit sequence is: ", DBHit)
            print("Hit seq in query starts from ", hit[2], " to ", hit[3], " and hit sequence is: ", QHit)
            print("Hit score: ", hit[4])
            hit_num += 1


def blast():
    read_database()
    remove_low_complexity_region()
    create_initial_seeds()
    make_permutation()
    find_Similar_seqs()
    print_Similar_Seqs()


w = int(input("Enter word_size: "))
word_Threshold = int(input("Enter word_Threshold: "))
HSP_Thershold = int(input("HSP_Thershold: "))
blast()
