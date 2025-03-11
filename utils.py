

_complement_trans_str = str.maketrans('acgtACGT', 'TGCATGCA')

genetic_code = {
  'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def dna_rev_comp(sequence: str) -> str:
  """
  Returns the uppercase reverse complement of the input DNA
  """
  return sequence.translate(_complement_trans_str)[::-1]


def translate_dna2aa(orf: str) -> str:
  """
  DNA to protein translation.
  Anything not from the genetic code is translated as "X".
  """
  protein = ''

  for i in range(0, (len(orf) // 3) * 3, 3):
    try:
      protein += genetic_code[orf[i:i + 3]]
    except KeyError:
      protein += 'X' #break
  return protein


def find_(string, 
         value_list
         ):
    """
    find the first occurence of a value in a list of values (e.g. read quality scores) in a string (e.g. read qualities)
    """
    
    indexes = [string.find(letter) for letter in value_list]
    #try: 
    ind = min([index for index in indexes if index != -1])
    # except:
    #     ind = 400
    return ind