#!/usr/bin/env python

import optparse
import sys
import os.path as op
import scipy.stats as ss

def hypergeom(m, n, n1, n2, p=False):
    """
    >>> hypergeom(1, 1000, 1000, 1000) # has to be shared.
    1.0

    >>> all(hypergeom(i, 1000, 1000, 1000) == 1.0 for i in range(100))
    True

    >>> hypergeom(1, 30000, 20, 20)
    0.013253396616299651

    >>> hypergeom(2, 30000, 20, 20)
    7.9649366037104485e-05

    >>> hypergeom(11, 30000, 20, 20)
    4.516176321800458e-11

    >>> hypergeom(10, 30000, 20, 20) # very low prob.
    4.516176321800458e-11

    >>> hypergeom(20, 30000, 20, 20) # very low chance that all are shared.
    4.516176321800458e-11

    """
    if m <= 0: return 1.0
    mmin = m - 1
    mmax = min(n1, n2)
    return ss.hypergeom.cdf(mmax, n, n1, n2) - ss.hypergeom.cdf(mmin, n, n1, n2)
    
def main():   

	file1=open("DHF6_significant_deltaP3_PPI201806_proteins.txt","r")
	file2=open("DiseaseModules.txt","r")
	file3=open("DHF6_significant_deltaP3_PPI201806_proteins_DiseaseModules.txt","w")
	
	disease=[]
	for line in file1:
		tmp=line.split()
		disease.append(tmp[0])
	file1.close()
	disease_set=set(disease)
	
	overlap=-1
	p=-1
	for line in file2:
		tmp=line.split('\t')
		endophenotype=[]
		file3.write("%s\t"%tmp[0])		
		for i in range(1,len(tmp)):
			endophenotype.append(tmp[i].strip())
		endophenotype_set=set(endophenotype)
		overlap=len(disease_set.intersection(endophenotype_set))
		p=hypergeom(overlap, 15489, len(disease),len(endophenotype))
		file3.write("%s\t%s\t%s\t"%(len(endophenotype),overlap,p))
		OG=list(disease_set.intersection(endophenotype_set))
		for t in range(overlap):
			file3.write("%s\t"%OG[t])
		file3.write("\n")
			
	file2.close()
	file3.close()
	
if __name__=='__main__':
	main()
		
    