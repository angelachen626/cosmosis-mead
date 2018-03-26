import os
import numpy as np

class Chainrun:
    
    # initialize the chain. input: inifile--cosmosis file to run, 
    # parfile--parameters of the chain, savefd--the folder to save all the results.
    # after instantiation, the object should have attributions inifile, parfile,
    # savefd, runfile(a new file cp from input file to avoid contaminaiton). runfile
    # will be in the same place as original ini file, with suffix '_chain' 
    def __init__(self,inifile='demo1', parfile=None, savefd='chaincosmosis'):
        self.inifile = inifile  #ini file
        self.parfile = parfile  #file of parameters to change in ini file
        self.savefd = savefd    #folder for results of every run in the chain
        self.runfile = self.inifile +'_chain.ini' #copy of ini file that we'll make change in
        os.system('cp '+self.inifile+'.ini '+self.runfile)
        opini = open(self.runfile,'r+')
        lines = opini.readlines()
        for line in lines:
            if line[:6]=='values':
                self.vlind = lines.index(line)
                self.vlfile = line[line.index('=')+1:].replace('\n','').replace(' ','')

            if line[:8] == 'save_dir':
                self.savind = lines.index(line)
        opini.close()
        os.system('cp '+self.vlfile+' '+self.vlfile[:-4]+'_chain'+'.ini ')
        self.vlfile =self.vlfile[:-4]+'_chain'+'.ini'
        self.chvlfile()
         
    
    # the method changing save_dir for each run.
    def chsavedir(self,savedir):
        opini = open(self.runfile,'r+')
        lines = opini.readlines()
        lines[self.savind] = 'save_dir='+savedir+'\n'
        print lines[self.savind]
        opini.seek(0)
        opini.truncate()
        opini.write(''.join(lines))
        opini.close()
        return 0

#changing value.ini file
    def chvlfile(self):
        opini = open(self.runfile,'r+')
        lines = opini.readlines()
        lines[self.vlind] = 'values='+self.vlfile+'\n'
        opini.seek(0)
        opini.truncate()
        opini.write(''.join(lines))
        opini.close()
        return 0

#changing parameters in values file
    def chvls(self,ind,value):
        opvl = open(self.vlfile,'r+')
        lines = opvl.readlines()
        chline = lines[ind]
        eqind = chline.index('=')
        chline = chline[:eqind+1]+str(value)+'\n'
        lines[ind] = chline
        opvl.seek(0)
        opvl.truncate()
        opvl.write(''.join(lines))
        opvl.close()
        return 0   

#give the object 1. the names [section,par name] of parameters to be changed. 2. the array of values of them. 3. the index of the line to be changed in the ini file.
    def findindofpar(self):
        oppar =  open(self.parfile,'r+')
        opvl = open(self.vlfile,'r+')
        parlines = oppar.readlines()
        vllines = opvl.readlines()
        oppar.close()
        opvl.close()
        self.parsnm= []
        self.parsvl=[]
        for line in parlines:
            if line[0] == '[':
                currentsection = line
            else:
                try:
                    eqind=line.index('=')
                except ValueError:
                    continue
                self.parsnm.append([currentsection,line[:eqind].replace(' ','')])#note currentsection has \n at the end, parname does not.
                self.parsvl.append(np.fromstring(line[eqind+1:],dtype=float,sep=' '))

        self.parsind= []

        for parnm in self.parsnm:
            conti = 0
            for line in vllines:
                if line[0]=='[':
                    currentsection = line
                    if conti:
                        conti=0
                if currentsection == parnm[0]:
                    conti = 1

                if line[:len(parnm[1])]==parnm[1] and conti:
                    self.parsind.append(vllines.index(line))
        
        
#    def chini():
        #replace the lines recorded in self.chdlines, a list
        
    
        
    def run_series(self):
        self.chsavedir(self.savefd)
        self.findindofpar()
        pargrid = np.array(np.meshgrid(*self.parsvl))
        parpts = zip(*np.reshape(pargrid,(len(self.parsvl),pargrid[0].size)))
        for parpt in parpts:
            savesuffix = '/chain'
            for i,value in enumerate(parpt):
                self.chvls(self.parsind[i],value)
                savesuffix = savesuffix+'_'+self.parsnm[i][0].replace('[','').replace(']','').replace('\n','')+'_'+self.parsnm[i][1]+'_'+str(value)
#                self.chsavedir(self.savefd+'aaa')
            self.chsavedir(self.savefd+savesuffix)
            os.system('cosmosis '+self.runfile)

    def cleanup(self):
        os.system('rm '+self.runfile)
        

def aaa():
    chain1 = Chainrun(inifile = 'demos/demo1',savefd = 'testdemo1',parfile='testdemo1/parfile.txt')
    chain1.run_series()
    print chain1.vlind, chain1.savind, str(chain1.parsind), str(chain1.parsnm), str(chain1.parsvl), chain1.runfile

if __name__ == "__main__":
    aaa()
