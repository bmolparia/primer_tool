"""
This module is a simple program that reads in a file with multiple sequences in fasta format and then
finds the primers that are required to join these sequences together using a yeast assembly protocol.
"""

import os, sys
import primer3

class InputError(Exception):
    pass

class Sequence(object):
    
    ## A sequence objectr to save all the attributes of the sequences

    def __init__(self,identifier, input_seq):
        
        self.seq  = input_seq
        self.iden = identifier

    def findRepeats(self):
        ## Function that looks for repeat regions within the sequence
        seq = self.seq
        pass

    def getTM(self):
        ## This function returns the melting temp of the sequence object
        seq = self.seq
        return primer3.calcTm(seq)
   
    @staticmethod
    def revComp(seq):
        ## Function to get the reverse complement of a sequence
        comp_bases = {'A':'T','G':'C','C':'G','T':'A'}
        new_seq = ''
        for base in seq:
            new_seq += comp_bases[base]
        return new_seq[::-1]
    
    
class parseInput(object):
    
    ## Handles the parsing of the input file 
        
    def __init__(self,input_file_path):

        self.inp = open(input_file_path,'U')

    def main(self):

        seqs = {}
        ## Initializing
        line = self.inp.readline()
        iden = ''
        tseq = 'A'
        cntr = 0
        
        while line:

            line = line.replace('\n','')
            if self.isNewSeq(line):
             
                try: 
                    if self.checkSeq(tseq):
                        seqO = Sequence(iden,tseq)  
                        seqs[cntr] = seqO
                except InputError,e:
                    print e
                    sys.exit()
                     
                ## re-initializing        
                temp = line.replace('>','').split('|')
                iden = temp[0]
                tseq = ''
                cntr += 1

            else:
                line = line.replace(' ','')
                tseq += line.upper()
            
            line = self.inp.readline()
       
        ## Adding in last read sequence
        seqO = Sequence(iden,tseq)
        seqs[cntr] = seqO
        
        seqs.pop(0) ## removing the 0th position sequence that was just used for initialization    
        
        return seqs
                               
    def isNewSeq(self,line):

        if line.startswith('>'):
            return True
        else:
            return False

    def checkSeq(self,seq):

        seq = list(set(list(seq)))
        for i in seq:
            if i not in ['A','T','G','C']:
                raise InputError('Sequence should just have A,T,G,C as base pairs. Found unrecognized character: '+i)

        return True


class Primers(object):

    def __init__(self,input_file_path, output_file_path, overlap_length=None,elongation_length=None):
        
        parseInst = parseInput(input_file_path)
        self.seqs = parseInst.main()            ## This is a dictionary containing all the sequences from the input file as sequence objects. 
                                                ## Every "key" is a number denoting the order of the sequences in the input
        self.outf = open(output_file_path,'w')

        if overlap_length == None:
            self.ovpL = 25
        else:
            self.ovpL = overlap_length

        if elongation_length == None:
            self.elonL = 55
        else:
            self.elonL = elongation_length
        
    def main(self):
        
        seq_dict = self.seqs
        num_seqs = max(seq_dict)
        ovp_len  = self.ovpL
        elon_len = self.elonL
        
        for ind in range(1,num_seqs+1):

        
            seqC  = seq_dict[ind] ## Current Sequence
            cseqC = Sequence.revComp(seqC.seq)
            
            primerF = ''
            primerR = ''

            if ind == 1:
                ## For first sequence
                seqN  = seq_dict[ind+1]                ## Next Sequence 
                cseqN = Sequence.revComp(seqN.seq)

                primerF = seqC.seq[:elon_len]
                primerR = cseqN[-ovp_len:] + cseqC[:elon_len]
                    
            elif ind == num_seqs:
                ## For last sequence
                seqP  = seq_dict[ind-1]                ## Previous Sequence

                primerF = seqP.seq[-ovp_len:] + seqC.seq[:elon_len]
                primerR = cseqC[:elon_len] 

            else:
                ## For all sequences in between
                seqP  = seq_dict[ind-1]                ## Previous Sequence
                seqN  = seq_dict[ind+1]                ## Next Sequence 
                cseqN = Sequence.revComp(seqN.seq)

                primerF = seqP.seq[-ovp_len:] + seqC.seq[:elon_len]
                primerR = cseqN[-ovp_len:] + cseqC[:elon_len] 

            Fobj = Sequence(seqC.iden,primerF)
            Robj = Sequence(seqC.iden,primerR)
            self.writeOut(Fobj,Robj)

            
        self.outf.close()

    def writeOut(self,Fobj,Robj):
        
        outp = self.outf
        FTm  = round(Fobj.getTM(),2)
        RTm  = round(Robj.getTM(),2)

        outp.write(Fobj.iden+'\n')
        outp.write(('\t').join( ['Forward Primer:',Fobj.seq, str(len(Fobj.seq)), str(FTm)] )+'\n')
        outp.write(('\t').join( ['Reverse Primer:',Robj.seq, str(len(Robj.seq)), str(RTm)] )+'\n\n')



if __name__ == '__main__':
   
    import argparse
    arg_parser = argparse.ArgumentParser(description="This module is a simple program that reads in a file with multiple sequences \
                                                      in fasta format and then finds the primers that are required to join these \
                                                      sequences together using a yeast assembly protocol.")
    
    arg_parser.add_argument('inputFile')
    arg_parser.add_argument('outputFile',default=os.path.join(os.getcwd(),'primers.txt'))
    arg_parser.add_argument('-oL','--overlap_length', type=int, default=None)
    arg_parser.add_argument('-eL','--elongation_length', type=int, default=None)
    
    args = arg_parser.parse_args()

    inp = args.inputFile
    out = args.outputFile
    oL  = args.overlap_length
    eL  = args.elongation_length
    
    pInst = Primers(inp,out,oL,eL)
    pInst.main()

