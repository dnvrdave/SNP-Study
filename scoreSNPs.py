# Python 2.7
# To print to a file (disables print command)
#from __future__ import print_function

# For Postgres database
import psycopg2

# Try/except block
import sys

# Debugger
import pdb

# Converts 'acgt' to '1234'
from string import maketrans

# Counts number of occurrences for each element of a list!
from collections import Counter

# To convert string to int
import ast

# Allows '==' operator to be passed to map f'n as 'eq'. Used in string comparison
#   Also needed to sort dict
import operator

# Profiler, for optimization
import cProfile

# regular expressions
import re

# To get a list of all files in a directory
from os import listdir
from os.path import isfile, join

def scoreSnps(trait):
  # Want False Alarm1 and False Alarm2 stats
  # e.g. hu0CF2EE	1	0 <-- has one alt allele,  did not report the trait == FalseAlarm1
  #      hu050E9C	2	0 <-- has two alt alleles, did not report the trait == FalseAlarm2
  # Connect to my local database, to get reported traits
  conme = psycopg2.connect("host='localhost' dbname='postgres' user='postgres' password='password'")
  curme = conme.cursor()
  print "Connected to my local database!\n"
  
  matchdir = 'C:/Users/David/Dropbox/BioLab/Disease-Gene Map/'
  indir    = 'C:/Users/David/Dropbox/BioLab/Disease-Gene Map/snpScores/' + trait + '/'
  # Get names of all files in the directory
  files = [f for f in listdir(indir) if isfile(join(indir, f))]
  rsid = []
  omatchFile = matchdir+trait+'_allScoresTrain.txt'
  ofile = open(omatchFile,'w')
  ofile.writelines('F1 has one alt allele,  did not report the trait\n')
  ofile.writelines('F2 has two alt alleles, did not report the trait\n')
  ofile.writelines('Miss has no alt alleles, but reported the trait\n')
  ofile.writelines('Hit1 has one alt allele,  and reported the trait\n')
  ofile.writelines('Hit2 has two alt alleles, and reported the trait\n\n')
  print ('F1 has one alt allele,  did not report the trait')
  print ('F2 has two alt alleles, did not report the trait')
  print ('Miss has no alt alleles, but reported the trait')
  print ('Hit1 has one alt allele,  and reported the trait')
  print ('Hit2 has two alt alleles, and reported the trait\n')
  print ('%14s\tBeta\tF1\tF2\tMiss\tHits\tHit1\tHit2\tHit0\tRep\tnotRep' % (trait))
  ofile.writelines('%14s\tBeta\tF1\tF2\tMiss\tHits\tHit1\tHit2\tHit0\t   Rep\tnotRep\n' % (trait))
  for i in range(len(files)):
    #print (files[i])
    rsidsplit = files[i].split("train.txt") # e.g. rs3019885train.txt
    rsid = rsidsplit[0]
    #SELECT beta FROM correl_gwas_train WHERE trait LIKE '%Asthma%' AND rsid = 'rs9319321';
    selectStr = "SELECT beta FROM correl_gwas_train WHERE trait LIKE '%" + trait + "%' AND rsid = '" + rsid + "';"
    #print ('%s' % (selectStr))

    try:
      curme.execute(selectStr)
    except:
      print "Can't SELECT\n"
    rows = curme.fetchone()
    if len(rows) > 0:
      beta = float(rows[0]) # rows is a tuple
      #print beta
      #pdb.set_trace()
    else:
      print "rsid not found\n"
      beta = 99
    fname = indir+files[i] # e.g. rs3019885train.txt
    f = open(fname, 'r')
    numReported = 0
    numNotReported = 0
    numHits = 0
    numHit0 = 0
    numHit1 = 0
    numHit2 = 0
    numMiss = 0
    numFalse1 = 0 # heterozygous
    numFalse2 = 0 # homozygous in alt allele
    for line in f:
        line = line.rstrip('\n')   # remove newline
        asplit = line.split("\t")  # parse line by tabs
        numAlt   = int(asplit[1])
        reported = int(asplit[2])
        if (reported):
          numReported+=1
          if (numAlt > 0):
            numHits+=1
            if (numAlt == 2):
              numHit2+=1
            else:
              numHit1+=1
          else:
            numMiss+=1
        else:    
          numNotReported+=1
          if (numAlt == 1):
            numFalse1+=1
          elif (numAlt == 2):
            numFalse2+=1
          else:
            numHit0+=1  # No Alt Alleles and Not Reported
    print ('%14s\t%4.2f\t%2s\t%2s\t%2s\t%2s\t%2s\t%2s\t%2s\t%3s\t%3s' % (rsid,beta,numFalse1,numFalse2,numMiss,numHits,numHit1,numHit2,numHit0,numReported,numNotReported))
    ofile.writelines('%14s\t%4.2f\t%2s\t%2s\t%4s\t%4s\t%4s\t%4s\t%4s\t%6s\t%6s\n' % (rsid,beta,numFalse1,numFalse2,numMiss,numHits,numHit1,numHit2,numHit0,numReported,numNotReported))    

def pull_allsnpsFast(trait):
  #
  if (0):
    files = ['hu011C57.vcf']
  else:
    files = ['hu011C57.vcf',  'hu016B28.vcf',           'hu0211D6.vcf',
           'hu025CEA.vcf',           'hu032C04.vcf',           'hu034DB1.vcf',
           'hu04DF3C.vcf',           'hu04F220.vcf',           'hu050E9C.vcf',
           'hu05FD49.vcf',           'hu089792.vcf',           'hu0CF2EE.vcf',
           'hu0D1FA1.vcf',           'hu0D879F.vcf',           'hu0E64A1.vcf',
           'hu1187FF.vcf',           'hu132B5C.vcf',           'hu1378E3.vcf',
           'hu15FECA.vcf',           'hu19C09F.vcf',           'hu1F73AB.vcf',
           'hu24C863.vcf',           'hu26B551.vcf',           'hu27FD1F.vcf',
           'hu2843C9.vcf',           'hu297562.vcf',           'hu2C1D94.vcf',
           'hu2DBF2D.vcf',           'hu2FEC01.vcf',           'hu3073E3.vcf',
           'hu33E2D9.vcf',           'hu34D5B9-GS01173.vcf',   'hu34D5B9-GS01670.vcf',
           'hu3C0611.vcf',           'hu3CAB43.vcf',           'hu3F864B.vcf',
           'hu4040B8.vcf',           'hu42D651.vcf',           'hu432EB5.vcf',
           'hu448C4B.vcf',           'hu44DCFF.vcf',           'hu470099.vcf',
           'hu48C4EB.vcf',           'hu4B0812.vcf',           'hu4BE6F2.vcf',
           'hu4BF398.vcf',           'hu4CA5B9.vcf',           'hu4FE0D1.vcf',
           'hu52B7E5.vcf',           'hu52F345.vcf',           'hu553620.vcf',
           'hu57A769.vcf',           'hu589D0B.vcf',           'hu599905.vcf',
           'hu5CD2C6.vcf',           'hu5E55F5.vcf',           'hu5FA322.vcf',
           'hu5FCE15.vcf',           'hu60180F.vcf',           'hu602487.vcf',
           'hu60AB7C.vcf',           'hu619F51.vcf',           'hu620F18.vcf',
           'hu627574.vcf',           'hu63EB0A.vcf',           'hu64DBF7.vcf',
           'hu661AD0.vcf',           'hu67EBB3.vcf',           'hu687B6B.vcf',
           'hu6C733E.vcf',           'hu6E4515.vcf',           'hu7123C1.vcf',
           'hu72C17A.vcf',           'hu76CAA5.vcf',           'hu775356.vcf',
           'hu7852C5.vcf',           'hu79F922.vcf',           'hu7A2F1D.vcf',
           'hu7A4AD1.vcf',           'hu7B594C.vcf',           'hu7DCBF9.vcf',
           'hu8229AE.vcf',           'hu82436A.vcf',           'hu82E689.vcf',
           'hu868880.vcf',           'hu8E2A35.vcf',           'hu8E87A9.vcf',
           'hu90B053.vcf',           'hu92C40A.vcf',           'hu92FD55.vcf',
           'hu939B7C.vcf',           'hu955EE1.vcf',           'huA05317.vcf',
           'huA0E089.vcf',           'huA49E22.vcf',           'huA4E2CF.vcf',
           'huA4F281.vcf',           'huAA245C.vcf',           'huAE4A11.vcf',
           'huAEADC0.vcf',           'huAEC1B0.vcf',           'huAFA81C.vcf',
           'huB1FD55.vcf',           'huB4883B.vcf',           'huB4940E.vcf',
           'huB4D223.vcf',           'huB4F9B2.vcf',           'huBA30D4.vcf',
           'huBAAC98.vcf',           'huBE0B25.vcf',           'huBEDA0B.vcf',
           'huC14AE1.vcf',           'huC29627.vcf',           'huC3160A.vcf', 'huE9B698.vcf']
  # Connect to my local database, to get reported traits
  conme = psycopg2.connect("host='localhost' dbname='postgres' user='postgres' password='password'")
  curme = conme.cursor()
  print "Connected to my local database!\n"
  
  # Use this to find huids that reported this trait
  queryProfiles = "SELECT id FROM reportedtraits WHERE trait LIKE '%"
  queryProfilesEnd = "%' ORDER BY id;"
  # Use this to get the rsids for the trait
  selectStr = "SELECT trait, chromosome, pos, tot_rsids_trait, rsid, pct_found_reported, tot_ids_reported, beta FROM correl_gwas_train WHERE trait LIKE '%"
  selectStrEnd = "%' ORDER BY chromosome::int, pos;"
  
  trait = trait.replace('\'','\'\'')  # Change ' to '' (e.g. Grave's --> Grave''s)
  traitWords = trait.split(" ")
  if len(traitWords) > 1:
        wordString = ""
        for i, tword in enumerate(traitWords):
          if i < len(traitWords)-1:
            # traitWords has more than one word, and this is not the last word
            # Don't confuse 'MIGRAINE WITHOUT AURA' as 'MIGRAINE WITH AURA'
            # Change this
            #   WHERE trait LIKE '%MIGRAINE%' AND trait LIKE '%
            # To this
            #   WHERE trait LIKE '%MIGRAINE WITH AURA%'
            if (tword.upper() == 'WITH') and ('MIGRAINE' in wordString):
              wordString = "MIGRAINE WITH AURA%' AND trait LIKE '%"
            else:
              #pdb.set_trace()
              wordString = wordString + tword.upper() + "%' AND trait LIKE '%"
          else:
            # the last word      
            wordString = wordString + tword.upper()
        #print wordString
        queryProfiles = queryProfiles + wordString + queryProfilesEnd
        #selectStr = selectStr + wordString + selectStrEnd
  else:
        # Single word trait      
        queryProfiles = queryProfiles + trait.upper() + queryProfilesEnd
  selectStr = selectStr + trait + selectStrEnd
  print queryProfiles
  print selectStr
  #pdb.set_trace()
  try:
        curme.execute(queryProfiles)
  except:
        print "Can't SELECT queryProfiles\n"
  # rowsp is list of all ids that reported the trait
  rowsp = curme.fetchall()
  totIdsReportedTrait = len(rowsp)
  print ('Total huids reporting trait %s %s\n' % (traitWords, totIdsReportedTrait))
  reportedDict = {}
  for i in range(len(files)):
    huidsplit = files[i].split(".")
    huid = huidsplit[0]
    huid = huid.replace('-GS01173','')
    huid = huid.replace('-GS01670','B')
    reportedDict[huid] = 0 # initialize
  for idd in rowsp:
    huidd = 'hu' + idd[0]
    print idd[0], huidd
    reportedDict[huidd] = 1   # key is huid, value is 1 if reported
  print ('%s' % (selectStr))
  try:
    curme.execute(selectStr)
  except:
    print "Can't SELECT\n"
  rows = curme.fetchall()
  if len(rows) > 0:
      traitDict    = {} # key=trait, value=trait number (alphabetical order)
      totrsidsDict = {} # key=trait, value=number of rsids found in GWAS for the trait
      pct_found_reportedDict = {}
      tot_reportedDict = {}
      traits    = []
      chromos   = []
      positions = []
      totrsids  = []
      rsids     = []
      pct_found_reported = []
      tot_reported = []
      beta_list = []
      chromopos = [] # x axis, one-up counter
      newchromo = [] # list of indexes where chromosome changes
      for i, row in enumerate(rows):
        # Skip if trait was not reported
        #if row[6] > 0:
        #if row[0].lower() == trait:
          traits.append(row[0])
          chromos.append(row[1])
          positions.append(row[2])
          totrsids.append(row[3])
          rsids.append(row[4])
          pct_found_reported.append(row[5])
          tot_reported.append(row[6])
          beta_list.append(row[7])
          chromopos.append(i)
  #pdb.set_trace()
  
  matchdir  = 'C:/Users/David/Dropbox/BioLab/Disease-Gene Map/snpScores/' + trait + '/'
  indir     = 'C:/Users/David/Data/PGP2015train/JustSNPs/'
  for j in range(len(rsids)):
    onesnp = rsids[j]
    omatchFile = matchdir+onesnp+'train.txt'
    ofile = open(omatchFile,'w')
    for i in range(len(files)):
      print (files[i])
      huidsplit = files[i].split(".")
      huid = huidsplit[0]
      huid = huid.replace('-GS01173','')
      huid = huid.replace('-GS01670','B')
      fname = indir+files[i]+'_allsnps.txt' # short snp lines from vcf file with both alleles
      #omatchFile = matchdir+files[i]+'_mthfr_snps.txt'
      f = open(fname, 'r')
      found  = 0
      numAlt = 0 # homozygous in ref allele (i.e. normal, no SNP)
      for line in f:
      #for j, line in enumerate(f):
        # rsid       refAllele   alt1   alt2
        #rs12455009	A	A	C
        line = line.rstrip('\n')   # remove newline
        asplit = line.split("\t")  # parse line by tabs
        #pdb.set_trace()
        if (asplit[0] == onesnp):
          # 1==1 Alt Allele
          # 2==2 Alt Alleles match
          # 3==2 Alt Alleles don't match (use 1 instead of 3 for auto analysis)
          # Assume sorted, i.e. if any alt allele matches the ref, it is the first one
          numAlt = 3 # both alleles are alternate and don't match each other
          if (asplit[1] == asplit[2]):
            numAlt = 1 # there is one alternate allele, heterozygous
          elif (asplit[2] == asplit[3]):
            numAlt = 2 # both alt alleles match, homozygous
          # huid RSID refAllele allele1 allele2
          #hu011C57	rs55780505  ATTTT       A	       ATTTT
          #pdb.set_trace()
          print ('found %16s\t%s\t%s\t%s\t%s\t%s\t%s' % (huid, asplit[0], asplit[1], asplit[2], asplit[3], numAlt, reportedDict[huid]))
          ofile.writelines('%s\t%s\t%s\n' % (huid, numAlt, reportedDict[huid]))
          found = 1
          break # end loop
      if (found == 0):    
        print ('not found %16s\t%s\t%s' % (huid, numAlt, reportedDict[huid]))
        ofile.writelines('%s\t%s\t%s\n' % (huid, numAlt, reportedDict[huid]))
    
def pull_onesnpFast(onesnp, trait):
  #
  if (0):
    files = ['hu011C57.vcf']
  else:
    files = ['hu011C57.vcf',  'hu016B28.vcf',           'hu0211D6.vcf',
           'hu025CEA.vcf',           'hu032C04.vcf',           'hu034DB1.vcf',
           'hu04DF3C.vcf',           'hu04F220.vcf',           'hu050E9C.vcf',
           'hu05FD49.vcf',           'hu089792.vcf',           'hu0CF2EE.vcf',
           'hu0D1FA1.vcf',           'hu0D879F.vcf',           'hu0E64A1.vcf',
           'hu1187FF.vcf',           'hu132B5C.vcf',           'hu1378E3.vcf',
           'hu15FECA.vcf',           'hu19C09F.vcf',           'hu1F73AB.vcf',
           'hu24C863.vcf',           'hu26B551.vcf',           'hu27FD1F.vcf',
           'hu2843C9.vcf',           'hu297562.vcf',           'hu2C1D94.vcf',
           'hu2DBF2D.vcf',           'hu2FEC01.vcf',           'hu3073E3.vcf',
           'hu33E2D9.vcf',           'hu34D5B9-GS01173.vcf',   'hu34D5B9-GS01670.vcf',
           'hu3C0611.vcf',           'hu3CAB43.vcf',           'hu3F864B.vcf',
           'hu4040B8.vcf',           'hu42D651.vcf',           'hu432EB5.vcf',
           'hu448C4B.vcf',           'hu44DCFF.vcf',           'hu470099.vcf',
           'hu48C4EB.vcf',           'hu4B0812.vcf',           'hu4BE6F2.vcf',
           'hu4BF398.vcf',           'hu4CA5B9.vcf',           'hu4FE0D1.vcf',
           'hu52B7E5.vcf',           'hu52F345.vcf',           'hu553620.vcf',
           'hu57A769.vcf',           'hu589D0B.vcf',           'hu599905.vcf',
           'hu5CD2C6.vcf',           'hu5E55F5.vcf',           'hu5FA322.vcf',
           'hu5FCE15.vcf',           'hu60180F.vcf',           'hu602487.vcf',
           'hu60AB7C.vcf',           'hu619F51.vcf',           'hu620F18.vcf',
           'hu627574.vcf',           'hu63EB0A.vcf',           'hu64DBF7.vcf',
           'hu661AD0.vcf',           'hu67EBB3.vcf',           'hu687B6B.vcf',
           'hu6C733E.vcf',           'hu6E4515.vcf',           'hu7123C1.vcf',
           'hu72C17A.vcf',           'hu76CAA5.vcf',           'hu775356.vcf',
           'hu7852C5.vcf',           'hu79F922.vcf',           'hu7A2F1D.vcf',
           'hu7A4AD1.vcf',           'hu7B594C.vcf',           'hu7DCBF9.vcf',
           'hu8229AE.vcf',           'hu82436A.vcf',           'hu82E689.vcf',
           'hu868880.vcf',           'hu8E2A35.vcf',           'hu8E87A9.vcf',
           'hu90B053.vcf',           'hu92C40A.vcf',           'hu92FD55.vcf',
           'hu939B7C.vcf',           'hu955EE1.vcf',           'huA05317.vcf',
           'huA0E089.vcf',           'huA49E22.vcf',           'huA4E2CF.vcf',
           'huA4F281.vcf',           'huAA245C.vcf',           'huAE4A11.vcf',
           'huAEADC0.vcf',           'huAEC1B0.vcf',           'huAFA81C.vcf',
           'huB1FD55.vcf',           'huB4883B.vcf',           'huB4940E.vcf',
           'huB4D223.vcf',           'huB4F9B2.vcf',           'huBA30D4.vcf',
           'huBAAC98.vcf',           'huBE0B25.vcf',           'huBEDA0B.vcf',
           'huC14AE1.vcf',           'huC29627.vcf',           'huC3160A.vcf', 'huE9B698.vcf']
  # Connect to my local database, to get reported traits
  conme = psycopg2.connect("host='localhost' dbname='postgres' user='postgres' password='password'")
  curme = conme.cursor()
  print "Connected to my local database!\n"
  
  # Use this to find huids that reported this trait
  queryProfiles = "SELECT id FROM reportedtraits WHERE trait LIKE '%"
  queryProfilesEnd = "%' ORDER BY id;"

  trait = trait.replace('\'','\'\'')  # Change ' to '' (e.g. Grave's --> Grave''s)
  traitWords = trait.split(" ")
  if len(traitWords) > 1:
        wordString = ""
        for i, tword in enumerate(traitWords):
          if i < len(traitWords)-1:
            # traitWords has more than one word, and this is not the last word
            # Don't confuse 'MIGRAINE WITHOUT AURA' as 'MIGRAINE WITH AURA'
            # Change this
            #   WHERE trait LIKE '%MIGRAINE%' AND trait LIKE '%
            # To this
            #   WHERE trait LIKE '%MIGRAINE WITH AURA%'
            if (tword.upper() == 'WITH') and ('MIGRAINE' in wordString):
              wordString = "MIGRAINE WITH AURA%' AND trait LIKE '%"
            else:
              #pdb.set_trace()
              wordString = wordString + tword.upper() + "%' AND trait LIKE '%"
          else:
            # the last word      
            wordString = wordString + tword.upper()
        #print wordString
        queryProfiles = queryProfiles + wordString + queryProfilesEnd
  else:
        # Single word trait      
        queryProfiles = queryProfiles + trait.upper() + queryProfilesEnd
  print queryProfiles
  #pdb.set_trace()
  try:
        curme.execute(queryProfiles)
  except:
        print "Can't SELECT queryProfiles\n"
  # rowsp is list of all ids that reported the trait
  rowsp = curme.fetchall()
  totIdsReportedTrait = len(rowsp)
  print ('Total huids reporting trait %s %s\n' % (traitWords, totIdsReportedTrait))
  reportedDict = {}
  for i in range(len(files)):
    huidsplit = files[i].split(".")
    huid = huidsplit[0]
    huid = huid.replace('-GS01173','')
    huid = huid.replace('-GS01670','B')
    reportedDict[huid] = 0 # initialize
  for idd in rowsp:
    huidd = 'hu' + idd[0]
    print idd[0], huidd
    reportedDict[huidd] = 1   # key is huid, value is 1 if reported
  #pdb.set_trace()
  
  matchdir  = 'C:/Users/David/Dropbox/BioLab/Disease-Gene Map/snpScores/' + trait + '/'
  indir     = 'C:/Users/David/Data/PGP2015train/JustSNPs/'
  rsid = []
  omatchFile = matchdir+onesnp+'train.txt'
  ofile = open(omatchFile,'w')
  for i in range(len(files)):
    print (files[i])
    huidsplit = files[i].split(".")
    huid = huidsplit[0]
    huid = huid.replace('-GS01173','')
    huid = huid.replace('-GS01670','B')
    fname = indir+files[i]+'_allsnps.txt' # short snp lines from vcf file with both alleles
    #omatchFile = matchdir+files[i]+'_mthfr_snps.txt'
    f = open(fname, 'r')
    found  = 0
    numAlt = 0 # homozygous in ref allele (i.e. normal, no SNP)
    for line in f:
    #for j, line in enumerate(f):
      # rsid       refAllele   alt1   alt2
      #rs12455009	A	A	C
      #if (line[0] != '#'): # Skip comment lines
      #if (line[0] != '#' and j<100000): # Skip comment lines. Stop early.
        line = line.rstrip('\n')   # remove newline
        asplit = line.split("\t")  # parse line by tabs
        #pdb.set_trace()
        if (asplit[0] == onesnp):
          # 1==1 Alt Allele
          # 2==2 Alt Alleles match
          # 3==2 Alt Alleles don't match (use 1 instead of 3 for auto analysis)
          # Assume sorted, i.e. if any alt allele matches the ref, it is the first one
          numAlt = 3 # both alleles are alternate and don't match each other
          if (asplit[1] == asplit[2]):
            numAlt = 1 # there is one alternate allele, heterozygous
          elif (asplit[2] == asplit[3]):
            numAlt = 2 # both alt alleles match, homozygous
          # huid RSID refAllele allele1 allele2
          #hu011C57	rs55780505  ATTTT       A	       ATTTT
          #pdb.set_trace()
          print ('found %16s\t%s\t%s\t%s\t%s\t%s\t%s' % (huid, asplit[0], asplit[1], asplit[2], asplit[3], numAlt, reportedDict[huid]))
          #ofile.writelines('%s\t%s\t%s\t%s\t%s\n' % (huid, asplit[0], asplit[1], asplit[2], asplit[3]))
          ofile.writelines('%s\t%s\t%s\n' % (huid, numAlt, reportedDict[huid]))
          found = 1
          break # end loop
    if (found == 0):    
      print ('not found %16s\t%s\t%s' % (huid, numAlt, reportedDict[huid]))
      ofile.writelines('%s\t%s\t%s\n' % (huid, numAlt, reportedDict[huid]))
    
# For each huid: numAltAlleles, reportedTheTrait, e.g. rs7216389train.txt, row1 = hu011C57	1	0
#pull_onesnpFast('rs7216389', 'asthma') # ASTHMA, chromo 17, Reads from allsnps files: rsid refAllele altAllele1 altAllele2
#pull_onesnpFast('rs3019885', 'asthma') # ASTHMA
#pull_allsnpsFast('Asthma') # Trait name must match correl_gwas_train case and format, e.g. Chrons_disease

# False alarm == 1 or 2 alt alleles but the person did not report the trait
#        Miss == 0 alt alleles but the person reported the trait
# Score: e.g. Asthma was reported by 30 people and had 48 SNPs in GWAS. Start with SNPs that had highest beta.
# For each SNP get a Miss score% and two False Alarm Scores and a Hit Score%:
# Miss Score% == 100x(numPeopleReportingWhoHadNoAltAllele/NumPeopleReporting)
# False Alarm Score1 == 100x(numPeopleWithOneAlleleWhoDidNotReport-numPeopleWithOneAlleleWhoDidReport)/numPeopleWithOneAllele
# False Alarm Score2 == 100x(numPeopleWithTwoAllelesWhoDidNotReport-numPeopleWithTwoAllelesWhoDidReport)/numPeopleWithTwoAlleles
# Hit Score% = 100x(numPeopleReportingWhoHad1or2AltAlleles + NumPeopleNotReportingWhoHadNoAltAllele)/totalPeople
# Read a set of SNP summary files like rs7216389train.txt and score each SNP
scoreSnps('Asthma') # Trait name must match correl_gwas_train case and format, e.g. Chrons_disease
