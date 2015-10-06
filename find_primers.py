"""
This is a simple program that reads in a file with multiple sequences in fasta format and then
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
                raise InputError('Sequence should just have A,T,G,C as bases. Found unrecognized character: '+i)

        return True


class Primers(object):

    def __init__(self,input_file_path, output_file_path, overlap_length=None,elongation_length=None):
        
        parseInst = parseInput(input_file_path)
        self.seqs = parseInst.main()            ## This is a dictionary containing all the sequences from the input file as sequence objects. 
                                                ## Every "key" is a number denoting the order of the sequences in the input
        self.outf = open(output_file_path,'w')
        self.junf = open(output_file_path+'.junction.txt','w')

        if overlap_length == None:
            self.ovpL = 25
        else:
            self.ovpL = overlap_length

        if elongation_length == None:
            self.elonL = 55
        else:
            self.elonL = elongation_length

        self.global_args ={'PRIMER_OPT_SIZE': 20,
                            'PRIMER_PICK_INTERNAL_OLIGO': 1,
                            'PRIMER_INTERNAL_MAX_SELF_END': 8,
                            'PRIMER_MIN_SIZE': 18,
                            'PRIMER_MAX_SIZE': 25,
                            'PRIMER_OPT_TM': 60.0,
                            'PRIMER_MIN_TM': 57.0,
                            'PRIMER_MAX_TM': 63.0,
                            'PRIMER_MIN_GC': 20.0,
                            'PRIMER_MAX_GC': 80.0,
                            'PRIMER_MAX_POLY_X': 100,
                            'PRIMER_INTERNAL_MAX_POLY_X': 100,
                            'PRIMER_SALT_MONOVALENT': 50.0,
                            'PRIMER_DNA_CONC': 50.0,
                            'PRIMER_MAX_NS_ACCEPTED': 0,
                            'PRIMER_MAX_SELF_ANY': 12,
                            'PRIMER_MAX_SELF_END': 8,
                            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                            'PRIMER_PAIR_MAX_COMPL_END': 8}
                    
    def main_primers(self):
        
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
    
    def junction_primers(self):
        
        seq_dict    = self.seqs
        num_seqs    = max(seq_dict)
        global_args = self.global_args

        outf = self.junf
        outf.write('## PrimerInfo = Primer_Seq\tTm\t%GC\n\n')
        ## Calculating the number of sets of multiplex primers and the number of primers in each set
        ## Example - if there are 11 sequences in total then the number of junction primers required = 10
        ## a total of 10 primers would do fine in this case. But if there are more, like 13, then a 6,6 split 
        ## would be better than a 10,2 split

        if (num_seqs-1)%10 != 0:
            segments = (num_seqs/10)+1
        else:
            segments = (num_seqs/10)
        
        set_length = (num_seqs-1)/segments
        remainder  = (num_seqs-1)%segments
        
        seg_arr = [set_length]*segments
        c = 1
        while c <= remainder:
            seg_arr[c-1] += 1
            c += 1
        #print seg_arr

        #Initializing the sequences and multiplex length 
        minimum_pcr_len = 40 # These are the minimum number of bp that will be a part of the PCR product from each segment
        main_counter  = 1
        setC = 1
        for seg in seg_arr:
            
            outf.write('###---------------------- SET '+str(setC)+'----------------------###\n')
            setC += 1

            seg_counter   = 0
            multiplex_len = 100
            while seg_counter < seg:

                seqC = seq_dict[main_counter].seq
                seqN = seq_dict[main_counter+1].seq
                
                temp_per_seg = multiplex_len-minimum_pcr_len
                pcr_template = seqC[-temp_per_seg:] + seqN[:temp_per_seg]
                
                global_args['PRIMER_PRODUCT_SIZE_RANGE'] = [[multiplex_len-10,multiplex_len+10]]
            
                seq_id   = seq_dict[main_counter].iden+' - '+seq_dict[main_counter+1].iden
                seq_args = {'SEQUENCE_ID':seq_id, 'SEQUENCE_TEMPLATE':pcr_template, 'SEQUENCE_INCLUDED_REGION': [0,len(pcr_template)]}
                
                res = primer3.bindings.designPrimers(seq_args, global_args)
        
                num_primer_pairs = res['PRIMER_PAIR_NUM_RETURNED']
            
                outf.write('\nJunction: '+seq_id+'\n\n')
                count_pair = 0
                while count_pair < num_primer_pairs:

                    left_primer   = str(res['PRIMER_LEFT_'+str(count_pair)+'_SEQUENCE'])
                    left_prim_tm  = str(res['PRIMER_LEFT_'+str(count_pair)+'_TM'])
                    left_prim_gc  = str(res['PRIMER_LEFT_'+str(count_pair)+'_GC_PERCENT'])
                    right_primer  = str(res['PRIMER_RIGHT_'+str(count_pair)+'_SEQUENCE'])
                    right_prim_tm = str(res['PRIMER_RIGHT_'+str(count_pair)+'_TM'])
                    right_prim_gc = str(res['PRIMER_RIGHT_'+str(count_pair)+'_GC_PERCENT'])
                    product_size  = str(res['PRIMER_PAIR_'+str(count_pair)+'_PRODUCT_SIZE'])
                    
                    outf.write( ('\t').join(['Forward: ', left_primer,  left_prim_tm,  left_prim_gc])  +'\n')
                    outf.write( ('\t').join(['Reverse: ', right_primer, right_prim_tm, right_prim_gc]) +'\n')
                    outf.write( 'Product Size: '+ product_size+'\n\n')
                    
                    count_pair += 1 
                
                seg_counter   += 1
                main_counter  += 1 
                multiplex_len += 50

        outf.close()

    def writeOut(self,Fobj,Robj):
        
        outp = self.outf
        FTm  = round(Fobj.getTM(),2)
        RTm  = round(Robj.getTM(),2)

        outp.write(Fobj.iden+'\n')
        outp.write(('\t').join( ['Forward Primer:',Fobj.seq, str(len(Fobj.seq)), str(FTm)] )+'\n')
        outp.write(('\t').join( ['Reverse Primer:',Robj.seq, str(len(Robj.seq)), str(RTm)] )+'\n\n')



if __name__ == '__main__':
   
    import argparse
    arg_parser = argparse.ArgumentParser(description="This is a simple program that reads in a file with multiple sequences \
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
    pInst.main_primers()
    pInst.junction_primers()

