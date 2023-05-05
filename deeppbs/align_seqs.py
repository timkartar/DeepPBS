import sys, os
from  Bio.SeqIO.PdbIO import PdbSeqresIterator
from Bio import pairwise2, SeqIO
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio import Align
from io import StringIO

import requests as r

def alignSeqs(seq1, seq2):
    aligner = Align.PairwiseAligner()

    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    alignments = aligner.align(seq1, seq2)
    alignment = list(alignments)[0]

    return str(alignment), alignment.score


if __name__ == "__main__":
    
    #dirname = "../../get/dataset/raw_pdbs_with_dna/"
    #for filename in os.listdir(dirname):
    #    for record in SeqIO.parse(os.path.join(dirname,filename),"pdb-seqres"):
    #        print(record.id, record.seq)
    
    dirname = "../../get/dataset/raw_pdbs_with_dna/" 
    for record in SeqIO.parse(os.path.join(dirname,"2r5y.pdb1"),"pdb-seqres"):
        if record.id == "A":
            break

    uniprot= "P09077" 
    baseUrl="http://www.uniprot.org/uniprot/"

    url = baseUrl+uniprot+".fasta"

    response = r.post(url)
    cData=''.join(response.text)

    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))[0].seq
    
    a, s = alignSeqs(record.seq,pSeq)

    print(a,s)


    
