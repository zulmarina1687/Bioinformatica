from Bio import SeqIO
import numpy as np

def datacleaning(path,dtype='fasta'):
    sequences = SeqIO.parse(path,dtype)
    for record in sequences:
        data = str(record.seq.upper())
    return data

class AligmentLocal:
    def __init__(self,_cadenaA,_cadenaB):
        self.matrixS = {'A':{'A':2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'C':{'A':-2,'C':2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'D':{'A':-2,'C':-2,'D':2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'E':{'A':-2,'C':-2,'D':-2,'E':2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'F':{'A':-2,'C':-2,'D':-2,'E':-2,'F':2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'G':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'H':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'I':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'K':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'L':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'M':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'N':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'P':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'Q':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'R':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'S':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':2,'T':-2,'V':-2,'W':-2,'Y':-2},
                        'T':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':2,'V':-2,'W':-2,'Y':-2},
                        'V':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':2,'W':-2,'Y':-2},
                        'W':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':2,'Y':-2},
                        'Y':{'A':-2,'C':-2,'D':-2,'E':-2,'F':-2,'G':-2,'H':-2,'I':-2,'K':-2,'L':-2,'M':-2,'N':-2,'P':-2,'Q':-2,'R':-2,'S':-2,'T':-2,'V':-2,'W':-2,'Y':2}}
        self.cadenaA = datacleaning(_cadenaA)
        self.cadenaB = datacleaning(_cadenaB)
        self.gap = 0
        self.matrixF = np.zeros((len(self.cadenaB)+1, len(self.cadenaA)+1))
    def updateMatrixS(self,match,missmatch):
        for key in self.matrixS:
            for key2 in self.matrixS[key]:
                if(key == key2):
                    self.matrixS[key][key2] = match
                else:
                    self.matrixS[key][key2] = missmatch
    def _local(self, _match=0, _missmatch=0, _gap=0):
        if(_match != 0 and _missmatch != 0):
            self.updateMatrixS(_match,_missmatch)
        self.gap = _gap
        self.makeFLocalMatrix()
        return self.localtraceback(_match, _missmatch, _gap)
    def localtraceback(self, match, missmatch, gap):
        vecA = []
        vecB = []
        newcadA = ""
        newcadB = ""
        _max_valor, _list_posiciones = self.makeFLocalMatrix()
        print(self.matrixF)
        num = len(_list_posiciones)/2 + 1
        num = int(num)
        for i in range (num):
            m=_list_posiciones[i][0]
            n=_list_posiciones[i][1]
            self.localBackTracking(m,n,newcadA,newcadB,vecA,vecB)
        print("Score:", self._score(match, missmatch, gap,vecA[0],vecB[0]))
        print("Alineamientos:")
        return vecA,vecB
    def makeFLocalMatrix(self):
        max_valor = 0
        list_posiciones = []
        filas, columnas = len(self.cadenaB), len(self.cadenaA)
        for i in range(filas):
            for j in range(columnas):
                self.matrixF[i+1][j+1] = max(
                    0, 
                    self.matrixF[i+1][j]+self.gap,
                    self.matrixF[i][j+1] + self.gap, 
                    self.matrixF[i][j] + self.matrixS[self.cadenaB[i]][self.cadenaA[j]])
                if self.matrixF[i+1][j+1] > max_valor:
                    max_valor = self.matrixF[i+1][j+1]
                    list_posiciones = []
                    list_posiciones.append([i,j])
                elif self.matrixF[i+1][j+1] == max_valor:
                    list_posiciones.append([i,j])
        return max_valor, list_posiciones
    def _score(self, match, missmatch, gap, cadA, cadB):
        score_=0
        for i in range(len(cadA)):
            if cadA[i] == cadB[i]:
                score_=score_ + match
            elif cadA[i]=='-' or cadB[i]=='-':
                score_=score_ + gap
            elif cadA[i] != cadB[i]:
                score_=score_ + missmatch
        return score_
    def localBackTracking(self,m,n,newcadA,newcadB,vecA,vecB):
        count = 0
        temporalA = ""
        temporalB = ""
        while( self.matrixF[m+1,n+1] ):
            score = self.matrixF[m+1,n+1]            
            scorediag = self.matrixF[m,n]
            scoreup = self.matrixF[m+1,n]
            scoreleft = self.matrixF[m,n+1]
            if score == scorediag + self.matrixS[self.cadenaB[m]][self.cadenaA[n]]:
                count = 1
                newcadA = self.cadenaA[n] + newcadA
                newcadB = self.cadenaB[m] + newcadB
                m-=1
                n-=1
            if scoreleft != 0 and score == scoreleft + self.gap:
                if count==1:
                    self.localBackTracking(m,n+1,"-"+temporalA,
                            self.cadenaB[m+1]+temporalB,vecA,vecB)
                else:
                    newcadB = self.cadenaB[m] + newcadB
                    newcadA = "-"+newcadA
                    m-=1
            if score == scoreup + self.gap:
                if count == 1:
                    self.localBackTracking(m+1,n,
                            self.cadenaA[n+1]+temporalA,
                            "-"+temporalB,vecA,vecB)
                else:
                    newcadB = "-" +newcadB
                    newcadA = self.cadenaA[n] + newcadA
                    n-=1
            count = 0
            temporalA = newcadA
            temporalB = newcadB
        vecA.append(newcadA)
        vecB.append(newcadB)