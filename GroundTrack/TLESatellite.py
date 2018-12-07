# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 16:13:05 2018

@author: ELias Obreque
"""

class Sat:
    def __init__(self, name, num, cla, YL, LNY, pl, EY, DEY, M1M,  M2M, BSTAR,
                 ET, esn, checksum, incl, RAAN, ec, ap, MA, nm, RNE):
        self.name   = name
        self.num    = num
        self.cla    = cla
        self.YL     = YL
        self.LNY    = LNY
        self.pl     = pl
        self.EY     = EY
        self.DEY    = DEY
        self.M1M    = M1M
        self.M2M    = M2M
        self.BSTAR  = BSTAR
        self.ET     = ET
        self.esn    = esn
        self.checksum = checksum
        self.incl   = incl
        self.RAAN   = RAAN
        self.ec     = ec
        self.ap     = ap
        self.MA     = MA
        self.nm     = nm
        self.RNE    = RNE
        print('Satellite created:', name)
              
    def ReadTLE(Archivo):
        TLE_open    = open(Archivo,'r')
        TLE_read    = TLE_open.read()
        i           = 0
        data        = ['','','']
        for line in TLE_read:
            if line == '\n':
                i += 1
            else:
                data[i] += line
        TLE_open.close()
        info        = [data[0],data[1],data[2]]
        return info 
    
    def setTLE(info):
        title = info[0]
        line_1 = info[1] 
        line_2 = info[2]
        satel  = Sat(title, line_1[2:7], line_1[8], line_1[9:11], line_1[11:14], line_1[14:17],
                     float('20'+line_1[18:20]), float(line_1[20:32]), line_1[33:43], line_1[44:52],
                     line_1[53:61], line_1[62], line_1[64:68], line_1[68], float(line_2[8:16]),
                     float(line_2[17:25]), float(line_2[26:33]), float(line_2[33:42]),
                     float(line_2[43:51]), float(line_2[52:63]), float(line_2[63:68]))
        return satel
    
    def getOE(self):
        n   = self.nm
        e   = self.ec/10000000.0
        i   = self.incl
        RAAN = self.RAAN
        w = self.ap
        M = self.MA
        return n, e, i, RAAN, w, M

        
        
        
    
